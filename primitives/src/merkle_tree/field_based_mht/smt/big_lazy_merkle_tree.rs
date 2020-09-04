use algebra::{ToBytes, to_bytes, FromBytes};

use crate::{
    crh::{FieldBasedHash, BatchFieldBasedHash, FieldBasedHashParameters},
    merkle_tree::{
        field_based_mht::{
            BatchFieldBasedMerkleTreeParameters,
            FieldBasedMerkleTreePath, FieldBasedMHTPath,
            smt::{
                Coord, OperationLeaf, ActionLeaf::Remove, BigMerkleTreeState, Error
            },
        },
    }
};

use rocksdb::{DB, Options};

use std::{
  collections::HashSet, marker::PhantomData, fs, path::Path,
};

#[derive(Debug)]
pub struct LazyBigMerkleTree<T: BatchFieldBasedMerkleTreeParameters> {
    // if unset, all DBs and tree internal state will be deleted when an instance of this struct
    // gets dropped
    persistent: bool,
    // path to state required to restore the tree
    state_path: Option<String>,
    // tree in-memory state
    state: BigMerkleTreeState<T>,
    // the height of the Merkle tree
    height: usize,
    // path to the db
    path_db: String,
    // stores the leaves
    database: DB,
    // path to the cache
    path_cache: String,
    // stores the cached nodes
    db_cache: DB,

    _parameters: PhantomData<T>,
}

impl<T: BatchFieldBasedMerkleTreeParameters> Drop for LazyBigMerkleTree<T> {
    fn drop(&mut self) {
        self.close()
    }
}

impl<T: BatchFieldBasedMerkleTreeParameters> LazyBigMerkleTree<T> {
    // Creates a new tree of specified `width`.
    // If `persistent` is specified, then DBs will be kept on disk and the tree state will be saved
    // so that the tree can be restored any moment later. Otherwise, no state will be saved on file
    // and the DBs will be deleted.
    pub fn new(
        width: usize,
        persistent: bool,
        state_path: Option<String>,
        path_db: String,
        path_cache: String
    ) -> Result<Self, Error> {
        let rate = <<<T::H as BatchFieldBasedHash>::BaseHash as FieldBasedHash>::Parameters as FieldBasedHashParameters>::R;
        assert_eq!(T::MERKLE_ARITY, 2); // For now we support only arity 2
        // Rate may also be smaller than the arity actually, but this assertion
        // is reasonable and simplify the design.
        assert_eq!(rate, T::MERKLE_ARITY);

        // If the tree must be persistent, than a path to which save the tree state must be
        // specified.
        if !persistent { assert!(state_path.is_none()) } else { assert!(state_path.is_some()) }

        let height = width as f64;
        let height = height.log(T::MERKLE_ARITY as f64) as usize;
        let state = BigMerkleTreeState::<T>::get_default_state(width, height);
        let path_db = path_db;
        let database = DB::open_default(path_db.clone())
            .map_err(|e| Error::Other(e.to_string()))?;
        let path_cache = path_cache;
        let db_cache = DB::open_default(path_cache.clone())
            .map_err(|e| Error::Other(e.to_string()))?;
        Ok(LazyBigMerkleTree {
            persistent,
            state_path,
            state,
            height,
            path_db,
            database,
            path_cache,
            db_cache,
            _parameters: PhantomData,
        })
    }

    // Restore a tree from a state saved at `state_path` and DBs at `path_db` and `path_cache`.
    // The new tree may be persistent or not, the actions taken in both cases are the same as
    // in `new` function.
    pub fn load(
        persistent: bool,
        state_path: String,
        path_db: String,
        path_cache: String
    ) -> Result<Self, Error> {
        let rate = <<<T::H as BatchFieldBasedHash>::BaseHash as FieldBasedHash>::Parameters as FieldBasedHashParameters>::R;
        assert_eq!(T::MERKLE_ARITY, 2); // For now we support only arity 2
        // Rate may also be smaller than the arity actually, but this assertion
        // is reasonable and simplify the design.
        assert_eq!(rate, T::MERKLE_ARITY);

        // Reads the state
        let state = { 
            let state_file = fs::File::open(state_path.clone())
                .map_err(|e| Error::Other(e.to_string()))?;
            BigMerkleTreeState::<T>::read(state_file)
        }.map_err(|e| Error::Other(e.to_string()))?;

        let height = state.width as f64;
        let height = height.log(T::MERKLE_ARITY as f64) as usize;

        let opening_options = Options::default();

        let path_db = path_db;
        let database = DB::open(&opening_options, path_db.clone())
            .map_err(|e| Error::Other(e.to_string()))?;

        let path_cache = path_cache;
        let db_cache = DB::open(&opening_options,path_cache.clone())
            .map_err(|e| Error::Other(e.to_string()))?;

        Ok(LazyBigMerkleTree {
            persistent,
            state_path: Some(state_path),
            state,
            height,
            path_db,
            database,
            path_cache,
            db_cache,
            _parameters: PhantomData,
        })
    }

    pub fn close(&mut self) {
        if !self.persistent {

            // Removes the state if present
            if self.state_path.is_some() && Path::new(&self.state_path.clone().unwrap()).exists() {
                match fs::remove_file(self.state_path.clone().unwrap()) {
                    Ok(_) => (),
                    Err(e) => println!("Error deleting tree state: {}", e)
                }
            }

            // Deletes the folder containing the db
            match fs::remove_dir_all(self.path_db.clone()) {
                Ok(_) => (),
                Err(e) => {
                    println!("Error deleting the folder containing the db: {}", e);
                }
            };
            // Deletes the folder containing the cache
            match fs::remove_dir_all(self.path_cache.clone()) {
                Ok(_) => (),
                Err(e) => {
                    println!("Error deleting the folder containing the db: {}", e);
                }
            };
        } else {
            // Don't delete the DBs and save required state on file in order to restore the
            // tree later.
            let tree_state_file = fs::File::create(self.state_path.clone().unwrap())
                .expect("Should be able to create file for tree state");

            self.state.write(tree_state_file)
                .expect("Should be able to write into tree state file");
        }
    }

    pub fn height(&self) -> usize { self.height }

    pub fn set_persistency(&mut self, persistency: bool) {
        self.persistent = persistency;
    }

    pub fn insert_to_cache(&self, coord: Coord, data:T::Data) {
        let elem = to_bytes!(data).unwrap();
        let index = bincode::serialize(&coord).unwrap();
        self.db_cache.put(index, elem).unwrap();
    }

    pub fn contains_key_in_cache(&self, coord:Coord) -> bool {
        let coordinates = bincode::serialize(&coord).unwrap();
        match self.db_cache.get(coordinates) {
            Ok(Some(_value)) => {
                return true;
            },
            Ok(None) => {
                return false;
            },
            Err(e) => {
                println!("operational problem encountered: {}", e);
                return false;
            },
        }
    }

    pub fn get_from_cache(&self, coord:Coord) -> Option<T::Data> {
        let coordinates = bincode::serialize(&coord).unwrap();
        match self.db_cache.get(coordinates) {
            Ok(Some(value)) => {
                let retrieved_elem = T::Data::read(value.as_slice()).unwrap();
                return Some(retrieved_elem);
            },
            Ok(None) => {
                return None;
            },
            Err(e) => {
                println!("operational problem encountered: {}", e);
                return None;
            },
        }
    }

    pub fn remove_from_cache(&self, coord: Coord) -> Option<T::Data>{
        let coordinates = bincode::serialize(&coord).unwrap();
        match self.db_cache.get(coordinates.clone()) {
            Ok(Some(value)) => {
                let retrieved_elem = T::Data::read(value.as_slice()).unwrap();
                let res = self.db_cache.delete(coordinates.clone());
                match res {
                    Ok(_) => {
                        return Some(retrieved_elem);
                    },
                    Err(e) => {
                        println!("Could not delete node from cache: {}", e);
                        return None;
                    }
                }
            },
            Ok(None) => {
                return None;
            },
            Err(e) => {
                println!("operational problem encountered: {}", e);
                return None;
            },
        }
    }

    pub fn insert_to_db(&self, idx: usize, data: T::Data) {
        let elem = to_bytes!(data).unwrap();
        let index = bincode::serialize(&idx).unwrap();
        self.database.put(index, elem).unwrap();
    }

    pub fn get_from_db(&self, idx: usize) -> Option<T::Data>{
        let index = bincode::serialize(&idx).unwrap();
        match self.database.get(index) {
            Ok(Some(value)) => {
                let retrieved_elem = T::Data::read(value.as_slice()).unwrap();
                return Some(retrieved_elem);
            },
            Ok(None) => {
                return None;
            },
            Err(e) => {
                println!("operational problem encountered: {}", e);
                return None;
            },
        }
    }

    pub fn remove_from_db(&self, idx: usize) -> Option<T::Data>{
        let index = bincode::serialize(&idx).unwrap();
        match self.database.get(index.clone()) {
            Ok(Some(value)) => {
                let retrieved_elem = T::Data::read(value.as_slice()).unwrap();
                let res = self.database.delete(index.clone());
                match res {
                    Ok(_) => {
                        return Some(retrieved_elem);
                    },
                    Err(e) => {
                        println!("Could not delete leaf from db: {}", e);
                        return None;
                    }
                }

            },
            Ok(None) => {
                return None;
            },
            Err(e) => {
                println!("operational problem encountered: {}", e);
                return None;
            },
        }
    }

    pub fn check_b_plus_caching_level_down(&mut self, coord: Coord) {

        let left_child_idx = coord.idx * T::MERKLE_ARITY;
        let left_child_height = 0;
        let left_coord = Coord { height: left_child_height, idx: left_child_idx };

        let right_child_idx = left_child_idx + 1;
        let right_child_height = 0;
        let right_coord = Coord { height: right_child_height, idx: right_child_idx };

        if (!self.state.present_node.contains(&left_coord)) | (!self.state.present_node.contains(&right_coord)) {
            if self.contains_key_in_cache(coord) {
                LazyBigMerkleTree::remove_node_from_cache(self, coord);
            }
        }
        return;
    }

    pub fn check_b_plus_caching(&mut self, coord: Coord) {
        assert_eq!(T::MERKLE_ARITY, 2, "Arity of the Merkle tree is not 2.");

        if coord.height <= 1 {
            return;
        }

        let left_child_idx = coord.idx * T::MERKLE_ARITY;
        let left_child_height = coord.height - 1;
        let left_coord = Coord { height: left_child_height, idx: left_child_idx };

        let right_child_idx = left_child_idx + 1;
        let right_child_height = left_child_height;
        let right_coord = Coord { height: right_child_height, idx: right_child_idx };

        if left_child_height == 1 {
            self.check_b_plus_caching_level_down(left_coord);
            self.check_b_plus_caching_level_down(right_coord);
            return;
        }
    }

    // Returns the hash value associated to the node.
    // If the node is in the cache, it retrieves from it.
    // If not, recomputes it.
    // Only used for nodes of level >= 1 (not for leaves).
    pub fn node(&mut self, coord: Coord) -> T::Data {

        assert_eq!(T::MERKLE_ARITY,2, "Arity of the Merkle tree is not 2.");

        //let coord = Coord{height, idx};
        // if the node is an empty node return the hash constant
        if !self.state.present_node.contains(&coord) {
            return T::EMPTY_HASH_CST[coord.height];
        }
        let res = self.get_from_cache(coord);

        // if not in the cache, recompute it
        if res == None {
            /* the node is not in the cache, compute it */
            let node_hash;
            if coord.height == 1 {
                /* get leaves to compute */
                let left_child_idx = coord.idx * T::MERKLE_ARITY;
                let left_child = self.get_from_db(left_child_idx);
                let left_hash: T::Data;
                if let Some(i) = left_child {
                    left_hash = i;
                } else {
                    left_hash = T::EMPTY_HASH_CST[0];
                }

                let right_child_idx = left_child_idx + 1;
                let right_child = self.get_from_db(right_child_idx);
                let right_hash: T::Data;
                if let Some(i) = right_child {
                    right_hash = i;
                } else {
                    right_hash = T::EMPTY_HASH_CST[0];
                }
                node_hash = Self::field_hash(&left_hash, &right_hash);
            } else {
                let height_child = coord.height - 1;
                let left_child_idx = coord.idx * T::MERKLE_ARITY;
                let coord_left = Coord { height: height_child, idx: left_child_idx };
                let left_child_hash = LazyBigMerkleTree::node(self, coord_left);

                let right_child_idx = left_child_idx + 1;
                let coord_right = Coord { height: height_child, idx: right_child_idx };
                let right_child_hash = LazyBigMerkleTree::node(self, coord_right);

                node_hash = Self::field_hash(&left_child_hash, &right_child_hash);
            }
            return node_hash;
        }

        res.unwrap()
    }

    pub fn get_root(&self) -> T::Data {
        self.state.root.clone()
    }

    pub fn remove_node_from_cache(&mut self, coord: Coord) {
        self.remove_from_cache(coord);
    }

    pub fn field_hash(x: &T::Data, y: &T::Data) -> T::Data{
        <<T::H as BatchFieldBasedHash>::BaseHash as FieldBasedHash>::init(None)
            .update(x.clone())
            .update(y.clone())
            .finalize()
    }

    pub fn batch_hash(input: &[T::Data]) -> Vec<T::Data> {
        <T::H as BatchFieldBasedHash>::batch_evaluate(input)
            .expect("Should be able to compute batch hash")
    }

    pub fn is_leaf_empty(&self, coord: Coord) -> bool {
        // check that the index of the leaf is less than the width of the Merkle tree
        assert!(coord.idx < self.state.width, "Leaf index out of bound.");
        // check that the coordinates of the node corresponds to the leaf level
        assert_eq!(coord.height, 0, "Coord of the node does not correspond to leaf level");

        !self.state.present_node.contains(&coord)
    }

    pub fn process_leaves (&mut self, vec_leaf_op: Vec<OperationLeaf<T::Data>>) -> T::Data {

        assert_eq!(T::MERKLE_ARITY, 2, "Arity of the Merkle tree is not 2.");

        // Validate inputs
        for i in 0..vec_leaf_op.len() {
            // check that the index of the leaf to be inserted/removed is in range
            assert!(vec_leaf_op[i].coord.idx < self.state.width, "Leaf index out of bound.");
        }

        // Mark nodes to recompute
        let mut visited_nodes:HashSet<Coord> = HashSet::new();
        let mut nodes_to_process_in_parallel:Vec<Vec<Coord>> = Vec::new();

        // process the first level
        let mut first_level_nodes:Vec<Coord> = Vec::new();
        for j in 0..vec_leaf_op.len() {
            let x = vec_leaf_op[j];
            let coord = x.coord;
            // visit parent
            let height_parent = coord.height + 1;
            let idx_parent = coord.idx / T::MERKLE_ARITY;
            let parent_coord = Coord {height: height_parent, idx: idx_parent};
            if !visited_nodes.contains(&parent_coord) {
                // parent node not visited yet
                visited_nodes.insert(parent_coord);
                // insert node to process in parallel
                first_level_nodes.push(parent_coord);
            }
        }
        nodes_to_process_in_parallel.push(first_level_nodes);

        // got to the upper levels until the root
        let mut height = 1;
        while height < self.height {

            visited_nodes.clear();
            let mut higher_level_nodes:Vec<Coord> = Vec::new();

            for j in 0..nodes_to_process_in_parallel[height - 1].len() {
                let coord = nodes_to_process_in_parallel[height - 1][j];
                let height_parent = coord.height + 1;
                let idx_parent = coord.idx / T::MERKLE_ARITY;
                let parent_coord = Coord { height: height_parent, idx: idx_parent };
                if !visited_nodes.contains(&parent_coord) {
                    // parent node not visited yet
                    visited_nodes.insert(parent_coord);
                    // insert node to process in parallel
                    higher_level_nodes.push(parent_coord);
                }
            }
            nodes_to_process_in_parallel.push(higher_level_nodes);
            height += 1;
        }

        // Inserts leaves in the db and marks as present/not present.
        for i in 0..vec_leaf_op.len() {
            let x = vec_leaf_op[i];
            let action = x.action;
            let coord = x.coord;
            let hash = x.hash;
            let idx = coord.idx;

            if action == Remove {
                self.remove_from_db(idx);
                self.state.present_node.remove(&coord);
            } else {
                self.insert_to_db(idx, hash.unwrap());
                self.state.present_node.insert(coord);
            }
        }

        // Compute hashes in parallel - first level
        let mut input_vec = Vec::new();
        let mut both_children_present = Vec::new();

        for j in 0..nodes_to_process_in_parallel[0].len() {
            let coord = nodes_to_process_in_parallel[0][j];

            let idx = coord.idx;
            let left_child_idx = idx * T::MERKLE_ARITY;
            let right_child_idx= left_child_idx + 1;
            let left_hash;
            let right_hash;
            let left_child_present;
            let right_child_present;
            let left_leaf = self.get_from_db(left_child_idx);
            if let Some(i) = left_leaf {
                left_hash = i;
                left_child_present = true;
            } else {
                left_hash = T::EMPTY_HASH_CST[0];
                left_child_present = false;
            }

            let right_leaf = self.get_from_db(right_child_idx);
            if let Some(i) = right_leaf {
                right_hash = i;
                right_child_present = true;
            } else {
                right_hash = T::EMPTY_HASH_CST[0];
                right_child_present = false;
            }
            input_vec.push(left_hash);
            input_vec.push(right_hash);
            if left_child_present && right_child_present {
                self.state.present_node.insert(coord);
            }
            if left_child_present && right_child_present {
                both_children_present.push(true);
            } else {
                both_children_present.push(false);
            }
        }
        // Process the input_vec using batch Poseidon hash
        let output_vec = Self::batch_hash(input_vec.as_slice());
        // Place the computed hash in a cache_parallel
        let mut index_output_vec = 0;
        for coord in nodes_to_process_in_parallel[0].clone() {
            self.state.cache_path.insert(coord, output_vec[index_output_vec]);
            if both_children_present[index_output_vec] == true {
                self.insert_to_cache(coord,output_vec[index_output_vec]);
            } else {
                self.remove_from_cache(coord);
            }
            index_output_vec += 1;
        }

        // Compute hashes in parallel - level > 1
        let mut height = 2;
        while height <= self.height {
            let mut input_vec = Vec::new();
            let mut both_children_present = Vec::new();
            for j in 0..nodes_to_process_in_parallel[height-1].len() {
                let coord = nodes_to_process_in_parallel[height -1][j];

                let idx = coord.idx;
                let left_child_idx = idx * T::MERKLE_ARITY;
                let right_child_idx = left_child_idx + 1;
                let left_hash: T::Data;
                let right_hash: T::Data;
                let left_child_coord = Coord { height: coord.height - 1, idx: left_child_idx};
                let right_child_coord = Coord { height: coord.height - 1, idx: right_child_idx};

                if self.state.cache_path.contains_key(&left_child_coord) {
                    left_hash = *self.state.cache_path.get(&left_child_coord).unwrap();
                } else {
                    left_hash = self.node(left_child_coord);
                }
                if self.state.cache_path.contains_key(&right_child_coord) {
                    right_hash = *self.state.cache_path.get(&right_child_coord).unwrap();
                } else {
                    right_hash = self.node(right_child_coord);
                }
                input_vec.push(left_hash);
                input_vec.push(right_hash);
                if self.state.present_node.contains(&left_child_coord) || self.state.present_node.contains(&right_child_coord){
                    self.state.present_node.insert(coord);
                } else {
                    self.state.present_node.remove(&coord);
                }
                if self.state.present_node.contains(&left_child_coord) && self.state.present_node.contains(&right_child_coord){
                    both_children_present.push(true);
                } else {
                    both_children_present.push(false);
                }
            }

            // Process the input_vec using batch Poseidon hash
            let output_vec = Self::batch_hash(input_vec.as_slice());

            // Place the computed hash in a cache_parallel
            let mut index_output_vec = 0;
            for coord in nodes_to_process_in_parallel[height-1].clone() {
                self.state.cache_path.insert(coord, output_vec[index_output_vec]);
                if both_children_present[index_output_vec] == true {
                    self.insert_to_cache(coord,output_vec[index_output_vec]);
                } else {
                    self.remove_from_cache(coord);
                    self.check_b_plus_caching(coord);
                }
                index_output_vec += 1;
            }

            height += 1;
        }

        self.state.root = *self.state.cache_path.get(&Coord{height:self.height,idx:0}).unwrap();
        self.state.root.clone()
    }

    // NB. Allows to get Merkle Path of empty leaves too
    pub fn get_merkle_path(&mut self, leaf_coord: Coord) -> FieldBasedMHTPath<<T::H as BatchFieldBasedHash>::BaseHash, T>
    {
        // check that the index of the leaf is less than the width of the Merkle tree
        assert!(leaf_coord.idx < self.state.width, "Leaf index out of bound.");

        // check that the coordinates of the node corresponds to the leaf level
        assert_eq!(leaf_coord.height, 0, "Coord of the node does not correspond to leaf level");

        let mut path = Vec::with_capacity(self.height);
        let mut node_idx = leaf_coord.idx;
        let mut height = 0;
        while height != self.height {

            // Estabilish if sibling is a left or right child
            let direction = node_idx % T::MERKLE_ARITY;
            let sibling_idx = if direction == 0 { node_idx + 1 } else { node_idx - 1 };

            // Get its hash
            let sibling_coord = Coord { height, idx: sibling_idx };
            let sibling =
                // If it's not empty
                if self.state.present_node.contains(&sibling_coord) {
                // If it's a leaf, it's in the DB and we can get it from there
                if height == 0 {
                    self.get_from_db(sibling_idx).unwrap()
                } else { // Otherwise, we need to recompute it
                    self.node(sibling_coord)
                }
            } else { // If it's empty then we can directly get the precomputed empty at this height
                T::EMPTY_HASH_CST[height]
            };

            // Push info to path
            path.push((vec![sibling], direction));

            height += 1; // go up one level
            node_idx = node_idx / T::MERKLE_ARITY; // compute the index of the parent
        }
        assert_eq!(self.node(Coord { height, idx: node_idx }), self.state.root);

        return FieldBasedMHTPath::<<T::H as BatchFieldBasedHash>::BaseHash, T>::new(path.as_slice());
    }
}

#[cfg(test)]
mod test {

    use algebra::{
        fields::{
          mnt4753::Fr as MNT4753Fr, mnt6753::Fr as MNT6753Fr
        },
        biginteger::BigInteger768,
        Field, UniformRand
    };

    use crate::{
        crh::{
            MNT4PoseidonHash, MNT6PoseidonHash
        },
        merkle_tree::field_based_mht::{
            naive::{
                NaiveFieldBasedMerkleTreeConfig, NaiveMerkleTree
            },
            smt::{
                OperationLeaf, Coord, ActionLeaf
            },
            poseidon::{
                MNT4PoseidonSMT, MNT4PoseidonSMTLazy,
                MNT6PoseidonSMT, MNT6PoseidonSMTLazy,
                MNT4753MHTPoseidonParameters, MNT6753MHTPoseidonParameters
            },
            FieldBasedMerkleTreeParameters, FieldBasedMerkleTreePath
        }
    };

    use std::{
        str::FromStr, path::Path, time::Instant
    };

    use rand::{
        rngs::OsRng, SeedableRng, RngCore
    };
    use rand_xorshift::XorShiftRng;

    struct MNT4753FieldBasedMerkleTreeParams;
    struct MNT6753FieldBasedMerkleTreeParams;

    impl NaiveFieldBasedMerkleTreeConfig for MNT4753FieldBasedMerkleTreeParams {
        const HEIGHT: usize = 6;
        type H = MNT4PoseidonHash;
    }

    impl NaiveFieldBasedMerkleTreeConfig for MNT6753FieldBasedMerkleTreeParams {
        const HEIGHT: usize = 6;
        type H = MNT6PoseidonHash;
    }

    type MNT4753FieldBasedMerkleTree = NaiveMerkleTree<MNT4753FieldBasedMerkleTreeParams>;
    type MNT6753FieldBasedMerkleTree = NaiveMerkleTree<MNT6753FieldBasedMerkleTreeParams>;

    #[test]
    fn process_leaves_mnt4_comp() {

        let num_leaves = 2usize.pow(23);
        let mut rng1 = XorShiftRng::seed_from_u64(9174123u64);
        let mut leaves_to_insert: Vec<OperationLeaf<MNT4753Fr>> = Vec::new();
        let mut leaves_to_remove: Vec<OperationLeaf<MNT4753Fr>> = Vec::new();

        let n = 1000;

        for _i in 0..n {
            let random: u64 = OsRng.next_u64();
            let idx = random % num_leaves as u64;
            let elem = MNT4753Fr::rand(&mut rng1);

            leaves_to_insert.push(OperationLeaf { coord: Coord { height: 0, idx: idx as usize }, action: ActionLeaf::Insert, hash: Some(elem.clone()) });
            leaves_to_remove.push(OperationLeaf { coord: Coord { height: 0, idx: idx as usize }, action: ActionLeaf::Remove, hash: None });
        }

        // Insertion

        let root1;
        let root2;
        let root3;
        let root4;
        {
            let mut smt1 = MNT4PoseidonSMT::new(
                num_leaves,
                false,
                None,
                String::from("./db_leaves_mnt4_comp_1"),
                String::from("./db_cache_mnt4_comp_1")
            ).unwrap();
            let leaves_to_process1 = leaves_to_insert.clone();
            let now = Instant::now();
            root1 = smt1.process_leaves_normal(leaves_to_process1);
            let new_now = Instant::now();

            let duration_normal = new_now.duration_since(now).as_millis();

            println!("duration normal = {}", duration_normal);

            // Removal

            let leaves_to_process3 = leaves_to_remove.clone();
            let now = Instant::now();
            root3 = smt1.process_leaves_normal(leaves_to_process3);
            let new_now = Instant::now();

            let duration_normal = new_now.duration_since(now).as_millis();

            println!("duration normal = {}",duration_normal);

        }

        {
            let mut smt2 = MNT4PoseidonSMTLazy::new(
                num_leaves,
                false,
                None,
                String::from("./db_leaves_mnt4_comp_2"),
                String::from("./db_cache_mnt4_comp_2")
            ).unwrap();
            let leaves_to_process2 = leaves_to_insert.clone();
            let now = Instant::now();
            root2 = smt2.process_leaves(leaves_to_process2);
            let new_now = Instant::now();

            let duration_fast = new_now.duration_since(now).as_millis();

            println!("duration fast = {}", duration_fast);

            let leaves_to_process4 = leaves_to_remove.clone();
            let now = Instant::now();
            root4 = smt2.process_leaves(leaves_to_process4);
            let new_now = Instant::now();

            let duration_fast = new_now.duration_since(now).as_millis();

            println!("duration fast = {}",duration_fast);

        }

        assert_eq!(root1, root2, "Roots are not equal");
        assert_eq!(root3, root4, "Roots are not equal");

        assert_eq!(root3, MNT4753MHTPoseidonParameters::EMPTY_HASH_CST[23], "Sequence of roots not equal");

    }

    #[test]
    fn process_leaves_mnt4() {

        let num_leaves = 32;
        let mut leaves_to_process: Vec<OperationLeaf<MNT4753Fr>> = Vec::new();

        leaves_to_process.push(OperationLeaf { coord: Coord { height: 0, idx: 0 }, action: ActionLeaf::Insert, hash: Some(MNT4753Fr::from_str("1").unwrap()) });
        leaves_to_process.push(OperationLeaf { coord: Coord { height: 0, idx: 9 }, action: ActionLeaf::Insert, hash: Some(MNT4753Fr::from_str("2").unwrap()) });
        leaves_to_process.push(OperationLeaf { coord: Coord { height: 0, idx: 16 }, action: ActionLeaf::Remove, hash: Some(MNT4753Fr::from_str("3").unwrap()) });
        leaves_to_process.push(OperationLeaf { coord: Coord { height: 0, idx: 29 }, action: ActionLeaf::Insert, hash: Some(MNT4753Fr::from_str("3").unwrap()) });
        leaves_to_process.push(OperationLeaf { coord: Coord { height: 0, idx: 16 }, action: ActionLeaf::Remove, hash: Some(MNT4753Fr::from_str("3").unwrap()) });

        let mut smt = MNT4PoseidonSMTLazy::new(
            num_leaves,
            false,
            None,
            String::from("./db_leaves_mnt4"),
            String::from("./db_cache_mnt4")
        ).unwrap();
        smt.process_leaves(leaves_to_process);

        //=============================================

        let mut leaves = Vec::new();
        leaves.push(MNT4753Fr::from_str("1").unwrap());
        for _ in 1..9 {
            let f = MNT4753Fr::zero();
            leaves.push(f);
        }
        leaves.push(MNT4753Fr::from_str("2").unwrap());
        for _ in 10..29 {
            let f = MNT4753Fr::zero();
            leaves.push(f);
        }
        leaves.push(MNT4753Fr::from_str("3").unwrap());
        for _ in 30..31 {
            let f = MNT4753Fr::zero();
            leaves.push(f);
        }
        let tree = MNT4753FieldBasedMerkleTree::new(&leaves).unwrap();

        assert_eq!(tree.root(), smt.state.root, "Roots are not equal");

    }

    #[test]
    fn process_leaves_mnt6() {

        let num_leaves = 2usize.pow(23);
        let mut rng1 = XorShiftRng::seed_from_u64(9174123u64);
        let mut leaves_to_insert: Vec<OperationLeaf<MNT6753Fr>> = Vec::new();
        let mut leaves_to_remove: Vec<OperationLeaf<MNT6753Fr>> = Vec::new();

        let n = 1000;

        for _i in 0..n {
            let random: u64 = OsRng.next_u64();
            let idx = random % num_leaves as u64;
            let elem = MNT6753Fr::rand(&mut rng1);

            leaves_to_insert.push(OperationLeaf { coord: Coord { height: 0, idx: idx as usize }, action: ActionLeaf::Insert, hash: Some(elem.clone()) });
            leaves_to_remove.push(OperationLeaf { coord: Coord { height: 0, idx: idx as usize }, action: ActionLeaf::Remove, hash: None });
        }

        // Insertion

        let root1;
        let root2;
        let root3;
        let root4;
        {
            let mut smt1 = MNT6PoseidonSMT::new(
                num_leaves,
                false,
                None,
                String::from("./db_leaves_mnt6_comp_1"),
                String::from("./db_cache_mnt6_comp_1")
            ).unwrap();
            let leaves_to_process1 = leaves_to_insert.clone();
            let now = Instant::now();
            root1 = smt1.process_leaves_normal(leaves_to_process1);
            let new_now = Instant::now();

            let duration_normal = new_now.duration_since(now).as_millis();

            println!("duration normal = {}", duration_normal);

            // Removal

            let leaves_to_process3 = leaves_to_remove.clone();
            let now = Instant::now();
            root3 = smt1.process_leaves_normal(leaves_to_process3);
            let new_now = Instant::now();

            let duration_normal = new_now.duration_since(now).as_millis();

            println!("duration normal = {}",duration_normal);

        }

        {
            let mut smt2 = MNT6PoseidonSMTLazy::new(
                num_leaves,
                false,
                None,
                String::from("./db_leaves_mnt6_comp_2"),
                String::from("./db_cache_mnt6_comp_2")
            ).unwrap();
            let leaves_to_process2 = leaves_to_insert.clone();
            let now = Instant::now();
            root2 = smt2.process_leaves(leaves_to_process2);
            let new_now = Instant::now();

            let duration_fast = new_now.duration_since(now).as_millis();

            println!("duration fast = {}", duration_fast);

            let leaves_to_process4 = leaves_to_remove.clone();
            let now = Instant::now();
            root4 = smt2.process_leaves(leaves_to_process4);
            let new_now = Instant::now();

            let duration_fast = new_now.duration_since(now).as_millis();

            println!("duration fast = {}",duration_fast);

        }

        assert_eq!(root1, root2, "Roots are not equal");
        assert_eq!(root3, root4, "Roots are not equal");

        assert_eq!(root3, MNT6753MHTPoseidonParameters::EMPTY_HASH_CST[23], "Sequence of roots not equal");

    }

    #[test]
    fn process_leaves_mnt6_comp() {

        let num_leaves = 32;
        let mut leaves_to_process: Vec<OperationLeaf<MNT6753Fr>> = Vec::new();

        leaves_to_process.push(OperationLeaf { coord: Coord { height: 0, idx: 0 }, action: ActionLeaf::Insert, hash: Some(MNT6753Fr::from_str("1").unwrap()) });
        leaves_to_process.push(OperationLeaf { coord: Coord { height: 0, idx: 9 }, action: ActionLeaf::Insert, hash: Some(MNT6753Fr::from_str("2").unwrap()) });
        leaves_to_process.push(OperationLeaf { coord: Coord { height: 0, idx: 16 }, action: ActionLeaf::Remove, hash: Some(MNT6753Fr::from_str("3").unwrap()) });
        leaves_to_process.push(OperationLeaf { coord: Coord { height: 0, idx: 29 }, action: ActionLeaf::Insert, hash: Some(MNT6753Fr::from_str("3").unwrap()) });
        leaves_to_process.push(OperationLeaf { coord: Coord { height: 0, idx: 16 }, action: ActionLeaf::Remove, hash: Some(MNT6753Fr::from_str("3").unwrap()) });

        let mut smt = MNT6PoseidonSMTLazy::new(
            num_leaves,
            false,
            None,
            String::from("./db_leaves_mnt6"),
            String::from("./db_cache_mnt6")).unwrap();
        smt.process_leaves(leaves_to_process);

        //=============================================

        let mut leaves = Vec::new();
        leaves.push(MNT6753Fr::from_str("1").unwrap());
        for _ in 1..9 {
            let f = MNT6753Fr::zero();
            leaves.push(f);
        }
        leaves.push(MNT6753Fr::from_str("2").unwrap());
        for _ in 10..29 {
            let f = MNT6753Fr::zero();
            leaves.push(f);
        }
        leaves.push(MNT6753Fr::from_str("3").unwrap());
        for _ in 30..31 {
            let f = MNT6753Fr::zero();
            leaves.push(f);
        }
        let tree = MNT6753FieldBasedMerkleTree::new(&leaves).unwrap();

        assert_eq!(tree.root(), smt.state.root, "Roots are not equal");
    }

    #[test]
    fn test_persistency() {
        let root = MNT4753Fr::new(
            BigInteger768([
                17131081159200801074,
                9006481350618111567,
                12051725085490156787,
                2023238364439588976,
                13194888104290656497,
                14162537977718443379,
                13575626123664189275,
                9267800406229717074,
                8973990559932404408,
                1830585533392189796,
                16667600459761825175,
                476991746583444
            ])
        );

        // create a persistent smt in a separate scope
        {
            let mut smt = MNT4PoseidonSMTLazy::new(
                32,
                true,
                Some(String::from("./persistency_test_info_lazy")),
                String::from("./db_leaves_persistency_test_info_lazy"),
                String::from("./db_cache_persistency_test_info_lazy")
            ).unwrap();

            //Insert some leaves in the tree
            let mut leaves_to_process: Vec<OperationLeaf<MNT4753Fr>> = Vec::new();

            leaves_to_process.push(OperationLeaf { coord: Coord { height: 0, idx: 0 }, action: ActionLeaf::Insert, hash: Some(MNT4753Fr::from_str("1").unwrap()) });
            leaves_to_process.push(OperationLeaf { coord: Coord { height: 0, idx: 9 }, action: ActionLeaf::Insert, hash: Some(MNT4753Fr::from_str("2").unwrap()) });

            smt.process_leaves(leaves_to_process);

            // smt gets dropped but its info should be saved
        }
        // files and directories should have been created
        assert!(Path::new("./persistency_test_info_lazy").exists());
        assert!(Path::new("./db_leaves_persistency_test_info_lazy").exists());
        assert!(Path::new("./db_cache_persistency_test_info_lazy").exists());

        // create a non-persistent smt in another scope by restoring the previous one
        {
            let mut smt = MNT4PoseidonSMTLazy::load(
                false,
                String::from("./persistency_test_info_lazy"),
                String::from("./db_leaves_persistency_test_info_lazy"),
                String::from("./db_cache_persistency_test_info_lazy")
            ).unwrap();

            // insert other leaves
            let mut leaves_to_process: Vec<OperationLeaf<MNT4753Fr>> = Vec::new();

            leaves_to_process.push(OperationLeaf { coord: Coord { height: 0, idx: 16 }, action: ActionLeaf::Insert, hash: Some(MNT4753Fr::from_str("10").unwrap()) });
            leaves_to_process.push(OperationLeaf { coord: Coord { height: 0, idx: 29 }, action: ActionLeaf::Insert, hash: Some(MNT4753Fr::from_str("3").unwrap()) });

            smt.process_leaves(leaves_to_process);

            // if truly state has been kept, then the equality below must pass, since `root` was
            // computed in one go with another smt
            assert_eq!(root, smt.state.root);

            // smt gets dropped and state and dbs are deleted
        }

        // files and directories should have been deleted
        assert!(!Path::new("./persistency_test_info_lazy").exists());
        assert!(!Path::new("./db_leaves_persistency_test_info_lazy").exists());
        assert!(!Path::new("./db_cache_persistency_test_info_lazy").exists());
    }

    #[test]
    fn merkle_tree_path_test_mnt4() {

        let num_leaves = 32;
        let mut leaves = Vec::with_capacity(num_leaves);
        let mut leaves_for_lazy_smt = Vec::with_capacity(num_leaves);
        let mut rng = XorShiftRng::seed_from_u64(1231275789u64);

        let mut smt = MNT4PoseidonSMTLazy::new(
            num_leaves,
            false,
            None,
            String::from("./db_leaves_merkle_tree_path_test_mnt4_lazy"),
            String::from("./db_cache_merkle_tree_path_test_mnt4_lazy")
        ).unwrap();

        // Generate random leaves, half of which empty
        for i in 0..num_leaves/2 {
            let leaf = MNT4753Fr::rand(&mut rng);
            leaves.push(leaf);
            leaves_for_lazy_smt.push(OperationLeaf { coord: Coord { height: 0, idx: i }, action: ActionLeaf::Insert, hash: Some(leaf)});
        }

        for i in num_leaves/2..num_leaves {
            let leaf = MNT4753Fr::zero();
            leaves.push(leaf);
            leaves_for_lazy_smt.push(OperationLeaf { coord: Coord { height: 0, idx: i }, action: ActionLeaf::Insert, hash: Some(leaf)});
        }

        // Compute the root of the tree, and do the same for a NaiveMHT, used here as reference
        let naive_tree = MNT4753FieldBasedMerkleTree::new(leaves.as_slice()).unwrap();
        let root = smt.process_leaves(leaves_for_lazy_smt);
        let naive_root = naive_tree.root();
        assert_eq!(root, naive_root);

        for i in 0..num_leaves {

            // Create and verify a FieldBasedMHTPath
            let path = smt.get_merkle_path(Coord{height: 0, idx: i});
            assert!(path.verify(6, &leaves[i], &root).unwrap());

            // Create and verify a Naive path
            let naive_path = naive_tree.generate_proof(i, &leaves[i]).unwrap();
            assert!(naive_path.verify(&naive_root, &leaves[i]).unwrap());

            // Assert the two paths are equal
            assert!(path.compare_with_binary(naive_path.path.as_slice()));
        }
    }

    #[test]
    fn merkle_tree_path_test_mnt6() {

        let num_leaves = 32;
        let mut leaves = Vec::with_capacity(num_leaves);
        let mut leaves_for_lazy_smt = Vec::with_capacity(num_leaves);
        let mut rng = XorShiftRng::seed_from_u64(1231275789u64);

        let mut smt = MNT6PoseidonSMTLazy::new(
            num_leaves,
            false,
            None,
            String::from("./db_leaves_merkle_tree_path_test_mnt6_lazy"),
            String::from("./db_cache_merkle_tree_path_test_mnt6_lazy")
        ).unwrap();

        // Generate random leaves, half of which empty
        for i in 0..num_leaves/2 {
            let leaf = MNT6753Fr::rand(&mut rng);
            leaves.push(leaf);
            leaves_for_lazy_smt.push(OperationLeaf { coord: Coord { height: 0, idx: i }, action: ActionLeaf::Insert, hash: Some(leaf)});
        }

        for i in num_leaves/2..num_leaves {
            let leaf = MNT6753Fr::zero();
            leaves.push(leaf);
            leaves_for_lazy_smt.push(OperationLeaf { coord: Coord { height: 0, idx: i }, action: ActionLeaf::Insert, hash: Some(leaf)});
        }

        // Compute the root of the tree, and do the same for a NaiveMHT, used here as reference
        let naive_tree = MNT6753FieldBasedMerkleTree::new(leaves.as_slice()).unwrap();
        let root = smt.process_leaves(leaves_for_lazy_smt);
        let naive_root = naive_tree.root();
        assert_eq!(root, naive_root);

        for i in 0..num_leaves {

            // Create and verify a FieldBasedMHTPath
            let path = smt.get_merkle_path(Coord{height: 0, idx: i});
            assert!(path.verify(6, &leaves[i], &root).unwrap());

            // Create and verify a Naive path
            let naive_path = naive_tree.generate_proof(i, &leaves[i]).unwrap();
            assert!(naive_path.verify(&naive_root, &leaves[i]).unwrap());

            // Assert the two paths are equal
            assert!(path.compare_with_binary(naive_path.path.as_slice()));
        }
    }
}