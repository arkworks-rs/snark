use algebra::{ToBytes, to_bytes, FromBytes};

use crate::{crh::{FieldBasedHash, BatchFieldBasedHash, FieldBasedHashParameters}, merkle_tree::{
    field_based_mht::{
        BatchFieldBasedMerkleTreeParameters, check_precomputed_parameters,
        FieldBasedMerkleTreePath, FieldBasedBinaryMHTPath,
        smt::{
            Coord, OperationLeaf, ActionLeaf::Remove, BigMerkleTreeState
        },
    },
}, ActionLeaf};

use rocksdb::{DB, Options, IteratorMode};

use std::{
  collections::HashSet, marker::PhantomData, fs, path::Path, io::{Error, ErrorKind}, sync::Arc,
};

#[derive(Debug)]
pub struct LazyBigMerkleTree<T: BatchFieldBasedMerkleTreeParameters> {
    // if unset, all DBs and tree internal state will be deleted when an instance of this struct
    // gets dropped
    pub(crate) persistent: bool,
    // tree in-memory state
    pub(crate) state: BigMerkleTreeState<T>,
    // the number of leaves
    pub(crate) width: usize,
    // path to the db
    pub(crate) path_db: String,
    // stores the leaves
    pub(crate) database: Arc<DB>,
    // path to the cache
    pub(crate) path_cache: String,
    // stores the cached nodes
    pub(crate) db_cache: Arc<DB>,

    _parameters: PhantomData<T>,
}

impl<T: BatchFieldBasedMerkleTreeParameters> Drop for LazyBigMerkleTree<T> {
    fn drop(&mut self) {
        self.flush()
    }
}

impl<T: BatchFieldBasedMerkleTreeParameters> LazyBigMerkleTree<T> {
    // Creates a new tree of specified `height`.
    // If `persistent` is specified, then DBs will be kept on disk and the tree state will be saved
    // so that the tree can be restored any moment later. Otherwise, no state will be saved on file
    // and the DBs will be deleted.
    pub fn new(
        height: usize,
        persistent: bool,
        path_db: String,
    ) -> Result<Self, Error> {
        assert!(check_precomputed_parameters::<T>(height));

        let rate = <<T::H  as FieldBasedHash>::Parameters as FieldBasedHashParameters>::R;

        assert_eq!(T::MERKLE_ARITY, 2); // For now we support only arity 2
        // Rate may also be smaller than the arity actually, but this assertion
        // is reasonable and simplify the design.
        assert_eq!(rate, T::MERKLE_ARITY);

        let state = BigMerkleTreeState::<T>::get_default_state(height);
        let width = T::MERKLE_ARITY.pow(height as u32);

        let path_db = path_db;
        let database = Arc::new(DB::open_default(path_db.clone())
            .map_err(|e| Error::new(ErrorKind::Other, e))?);

        let mut path_cache = path_db.clone();
        path_cache.push_str("_cache");
        let db_cache = Arc::new(DB::open_default(path_cache.clone())
            .map_err(|e| Error::new(ErrorKind::Other, e))?);

        Ok(Self {
            persistent,
            state,
            width,
            path_db,
            database,
            path_cache,
            db_cache,
            _parameters: PhantomData,
        })
    }

    // Restore a tree from a state saved at `state_path` and DBs at `path_db` and `path_cache`.
    // The new tree may be persistent or not, the actions taken in both cases are the same as
    // in `new` function. The function is able to restore the tree iff at least `path_db` is
    // present.
    pub fn load(
        height: usize,
        persistent: bool,
        path_db: String,
    ) -> Result<Self, Error> {

        let leaves_db_path = Path::new(path_db.as_str());
        let cache_db_str = {
            let mut t = path_db.clone();
            t.push_str("_cache");
            t
        };
        let cache_db_path = Path::new(cache_db_str.as_str());

        // Temporary: clean up anyway, as we have not implemented still a way
        // to check for dbs consistency
        if cache_db_path.exists() { fs::remove_dir_all(cache_db_path)? }

        // The DB containing all the leaves it's enough to reconstruct a consistent tree.
        // Even if other info are available on disk, we would still need to check that
        // they stand for the same tree state, which is expensive nevertheless (and still
        // requires to compute the root). At this point, as a temporary solution, it's worth
        // in any case to reconstruct the whole tree given the leaves.
        if leaves_db_path.exists()
        {
            Self::restore_from_leaves(height, persistent, path_db)
        } else
        {
            // If `database` has been lost, it's not possible to load or recover anything
            Err(Error::new(
                ErrorKind::InvalidData,
                "Unable to restore MerkleTree: leaves not available")
            )
        }
    }

    // Use the database at db_path to restore the tree, thus recreating `state` and `db_cache`.
    // This operation is expensive.
    fn restore_from_leaves(
        height: usize,
        persistent: bool,
        path_db: String,
    ) -> Result<Self, Error> {

        let rate = <<T::H as FieldBasedHash>::Parameters as FieldBasedHashParameters>::R;
        assert_eq!(T::MERKLE_ARITY, 2); // For now we support only arity 2
        // Rate may also be smaller than the arity actually, but this assertion
        // is reasonable and simplify the design.
        assert_eq!(rate, T::MERKLE_ARITY);

        // Create new state
        let state = BigMerkleTreeState::<T>::get_default_state(height);
        assert!(check_precomputed_parameters::<T>(state.height));
        let width = T::MERKLE_ARITY.pow(height as u32);

        let opening_options = Options::default();

        // Restore leaves DB
        let path_db = path_db.clone();
        let database = Arc::new(DB::open(&opening_options, path_db.clone())
            .map_err(|e| Error::new(ErrorKind::Other, e))?);
        // It's fine here to have one instance writing and the other reading: first because
        // they are not concurrent and second because `database` already contain all the leaves
        // we are going to add, so there will be no writings at all.
        let database_read = database.clone();

        // Create new cache DB
        let mut path_cache = path_db.clone();
        path_cache.push_str("_cache");
        let db_cache = Arc::new(DB::open_default(path_cache.clone())
            .map_err(|e| Error::new(ErrorKind::Other, e))?);

        // Create new tree instance
        let mut new_tree = Self {
            persistent,
            state,
            width,
            path_db,
            database,
            path_cache,
            db_cache,
            _parameters: PhantomData,
        };

        // Re-add the leaves to the tree in order to update state and db_cache
        let db_iter = database_read.iterator(IteratorMode::Start);
        let mut leaves_to_insert = Vec::with_capacity(db_iter.size_hint().0);

        // For the moment, let's load all the leaves in memory and update the tree
        for (leaf_index, leaf_val) in db_iter {
            let index = u32::read(&*leaf_index)?;
            let leaf = T::Data::read(&*leaf_val)?;
            leaves_to_insert.push(OperationLeaf {
                coord: Coord { height: 0, idx: index as usize },
                action: ActionLeaf::Insert,
                hash: Some(leaf.clone())
            });
        }

        new_tree.process_leaves(leaves_to_insert.as_slice());
        Ok(new_tree)
    }

    pub fn flush(&mut self) {

        if !self.persistent {

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
        }
    }
    
    pub fn set_persistency(&mut self, persistency: bool) {
        self.persistent = persistency;
    }

    fn insert_to_cache(&self, coord: Coord, data:T::Data) {
        let elem = to_bytes!(data).unwrap();
        let index = to_bytes!(coord).unwrap();
        self.db_cache.put(index, elem).unwrap();
    }

    fn contains_key_in_cache(&self, coord:Coord) -> bool {
        let coordinates = to_bytes!(coord).unwrap();
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

    fn get_from_cache(&self, coord:Coord) -> Option<T::Data> {
        let coordinates = to_bytes!(coord).unwrap();
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

    fn remove_from_cache(&self, coord: Coord) -> Option<T::Data>{
        let coordinates = to_bytes!(coord).unwrap();
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

    fn insert_to_db(&self, idx: usize, data: T::Data) {
        let elem = to_bytes!(data).unwrap();
        let index = to_bytes!(idx as u32).unwrap();
        self.database.put(index, elem).unwrap();
    }

    fn get_from_db(&self, idx: usize) -> Option<T::Data>{
        let index = to_bytes!(idx as u32).unwrap();
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

    fn remove_from_db(&self, idx: usize) -> Option<T::Data>{
        let index = to_bytes!(idx as u32).unwrap();
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

    fn check_b_plus_caching_level_down(&mut self, coord: Coord) {

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

    fn check_b_plus_caching(&mut self, coord: Coord) {
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
    fn node(&mut self, coord: Coord) -> T::Data {

        assert_eq!(T::MERKLE_ARITY,2, "Arity of the Merkle tree is not 2.");

        //let coord = Coord{height, idx};
        // if the node is an empty node return the hash constant
        if !self.state.present_node.contains(&coord) {
            return T::EMPTY_HASH_CST.unwrap().nodes[coord.height];
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
                    left_hash = T::EMPTY_HASH_CST.unwrap().nodes[0];
                }

                let right_child_idx = left_child_idx + 1;
                let right_child = self.get_from_db(right_child_idx);
                let right_hash: T::Data;
                if let Some(i) = right_child {
                    right_hash = i;
                } else {
                    right_hash = T::EMPTY_HASH_CST.unwrap().nodes[0];
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

    pub fn height(&self) -> usize { self.state.height }

    fn remove_node_from_cache(&mut self, coord: Coord) {
        self.remove_from_cache(coord);
    }

    fn field_hash(x: &T::Data, y: &T::Data) -> T::Data{
        <T::H as FieldBasedHash>::init(None)
            .update(x.clone())
            .update(y.clone())
            .finalize()
    }

    fn batch_hash(input: &[T::Data]) -> Vec<T::Data> {
        <T::BH as BatchFieldBasedHash>::batch_evaluate(input)
            .expect("Should be able to compute batch hash")
    }

    pub fn is_leaf_empty(&self, coord: Coord) -> bool {
        // check that the index of the leaf is less than the width of the Merkle tree
        assert!(coord.idx < self.width, "Leaf index out of bound.");
        // check that the coordinates of the node corresponds to the leaf level
        assert_eq!(coord.height, 0, "Coord of the node does not correspond to leaf level");

        !self.state.present_node.contains(&coord)
    }

    pub fn process_leaves (&mut self, vec_leaf_op: &[OperationLeaf<T::Data>]) -> T::Data {

        assert_eq!(T::MERKLE_ARITY, 2, "Arity of the Merkle tree is not 2.");

        // Validate inputs
        for i in 0..vec_leaf_op.len() {
            // check that the index of the leaf to be inserted/removed is in range
            assert!(vec_leaf_op[i].coord.idx < self.width, "Leaf index out of bound.");
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
        while height < self.state.height {

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

            let mut left_child_present = true;
            let left_hash = self.get_from_db(left_child_idx).unwrap_or_else(|| {
                left_child_present = false;
                T::EMPTY_HASH_CST.unwrap().nodes[0]
            });

            let mut right_child_present = true;
            let right_hash = self.get_from_db(right_child_idx).unwrap_or_else(|| {
                right_child_present = false;
                T::EMPTY_HASH_CST.unwrap().nodes[0]
            });

            input_vec.push(left_hash);
            input_vec.push(right_hash);

            if left_child_present || right_child_present {
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
            if both_children_present[index_output_vec] {
                self.insert_to_cache(coord,output_vec[index_output_vec]);
            } else {
                self.remove_from_cache(coord);
            }
            index_output_vec += 1;
        }

        // Compute hashes in parallel - level > 1
        let mut height = 2;
        while height <= self.state.height {
            let mut input_vec = Vec::new();
            let mut both_children_present = Vec::new();
            for j in 0..nodes_to_process_in_parallel[height-1].len() {
                let coord = nodes_to_process_in_parallel[height -1][j];

                let idx = coord.idx;
                let left_child_idx = idx * T::MERKLE_ARITY;
                let right_child_idx = left_child_idx + 1;

                let left_child_coord = Coord { height: coord.height - 1, idx: left_child_idx};
                let right_child_coord = Coord { height: coord.height - 1, idx: right_child_idx};

                let left_hash = if self.state.cache_path.contains_key(&left_child_coord) {
                    *self.state.cache_path.get(&left_child_coord).unwrap()
                } else {
                    self.node(left_child_coord)
                };

                let right_hash = if self.state.cache_path.contains_key(&right_child_coord) {
                    *self.state.cache_path.get(&right_child_coord).unwrap()
                } else {
                    self.node(right_child_coord)
                };
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

        self.state.root = *self.state.cache_path.get(&Coord{height:self.state.height,idx:0}).unwrap();
        self.state.root.clone()
    }

    // NB. Allows to get Merkle Path of empty leaves too
    pub fn get_merkle_path(&mut self, leaf_coord: Coord) -> FieldBasedBinaryMHTPath<T>
    {
        // check that the index of the leaf is less than the width of the Merkle tree
        assert!(leaf_coord.idx < self.width, "Leaf index out of bound.");

        // check that the coordinates of the node corresponds to the leaf level
        assert_eq!(leaf_coord.height, 0, "Coord of the node does not correspond to leaf level");

        let mut path = Vec::with_capacity(self.state.height);
        let mut node_idx = leaf_coord.idx;
        let mut height = 0;
        while height != self.state.height {

            // Estabilish if sibling is a left or right child
            let (sibling_idx, direction) = if node_idx % T::MERKLE_ARITY == 0 {
                (node_idx + 1, false)
            } else {
                (node_idx - 1, true)
            };

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
                T::EMPTY_HASH_CST.unwrap().nodes[height]
            };

            // Push info to path
            path.push((sibling, direction));

            height += 1; // go up one level
            node_idx = node_idx / T::MERKLE_ARITY; // compute the index of the parent
        }
        assert_eq!(self.node(Coord { height, idx: node_idx }), self.state.root);

        return FieldBasedBinaryMHTPath::<T>::new(path);
    }
}

#[cfg(test)]
mod test {

    use algebra::{
        fields::{
          mnt4753::Fr as MNT4753Fr, mnt6753::Fr as MNT6753Fr
        },
        biginteger::BigInteger768,
        Field, UniformRand,
        ToBytes, to_bytes, FromBytes,
    };

    use crate::{
        crh::parameters::{MNT4PoseidonHash, MNT4BatchPoseidonHash, MNT6PoseidonHash, MNT6BatchPoseidonHash},
        merkle_tree::field_based_mht::{
            naive:: NaiveMerkleTree,
            smt::{OperationLeaf, Coord, ActionLeaf, BigMerkleTree, LazyBigMerkleTree},
            parameters::{
                MNT4753_MHT_POSEIDON_PARAMETERS, MNT6753_MHT_POSEIDON_PARAMETERS
            },
            FieldBasedMerkleTreeParameters, BatchFieldBasedMerkleTreeParameters,
            FieldBasedMerkleTreePrecomputedEmptyConstants, FieldBasedMerkleTreePath,
            FieldBasedBinaryMHTPath,
        },

    };

    use std::{
        str::FromStr, path::Path, time::Instant
    };

    use rand::{
        rngs::OsRng, SeedableRng, RngCore
    };
    use rand_xorshift::XorShiftRng;

    #[derive(Clone, Debug)]
    struct MNT4753FieldBasedMerkleTreeParams;
    impl FieldBasedMerkleTreeParameters for MNT4753FieldBasedMerkleTreeParams {
        type Data = MNT4753Fr;
        type H = MNT4PoseidonHash;
        const MERKLE_ARITY: usize = 2;
        const EMPTY_HASH_CST: Option<FieldBasedMerkleTreePrecomputedEmptyConstants<'static, Self::H>> = Some(MNT4753_MHT_POSEIDON_PARAMETERS);
    }
    impl BatchFieldBasedMerkleTreeParameters for MNT4753FieldBasedMerkleTreeParams {
        type BH = MNT4BatchPoseidonHash;
    }
    type MNT4753FieldBasedMerkleTree = NaiveMerkleTree<MNT4753FieldBasedMerkleTreeParams>;
    type MNT4PoseidonSMT = BigMerkleTree<MNT4753FieldBasedMerkleTreeParams>;
    type MNT4PoseidonSMTLazy = LazyBigMerkleTree<MNT4753FieldBasedMerkleTreeParams>;
    type MNT4MerklePath = FieldBasedBinaryMHTPath<MNT4753FieldBasedMerkleTreeParams>;

    #[derive(Clone, Debug)]
    struct MNT6753FieldBasedMerkleTreeParams;
    impl FieldBasedMerkleTreeParameters for MNT6753FieldBasedMerkleTreeParams {
        type Data = MNT6753Fr;
        type H = MNT6PoseidonHash;
        const MERKLE_ARITY: usize = 2;
        const EMPTY_HASH_CST: Option<FieldBasedMerkleTreePrecomputedEmptyConstants<'static, Self::H>> = Some(MNT6753_MHT_POSEIDON_PARAMETERS);
    }
    impl BatchFieldBasedMerkleTreeParameters for MNT6753FieldBasedMerkleTreeParams {
        type BH = MNT6BatchPoseidonHash;
    }
    type MNT6753FieldBasedMerkleTree = NaiveMerkleTree<MNT6753FieldBasedMerkleTreeParams>;
    type MNT6PoseidonSMT = BigMerkleTree<MNT6753FieldBasedMerkleTreeParams>;
    type MNT6PoseidonSMTLazy = LazyBigMerkleTree<MNT6753FieldBasedMerkleTreeParams>;
    type MNT6MerklePath = FieldBasedBinaryMHTPath<MNT6753FieldBasedMerkleTreeParams>;

    const TEST_HEIGHT_1: usize = 5;

    #[test]
    fn process_leaves_mnt4() {

        let mut leaves_to_process: Vec<OperationLeaf<MNT4753Fr>> = Vec::new();

        leaves_to_process.push(OperationLeaf { coord: Coord { height: 0, idx: 0 }, action: ActionLeaf::Insert, hash: Some(MNT4753Fr::from_str("1").unwrap()) });
        leaves_to_process.push(OperationLeaf { coord: Coord { height: 0, idx: 9 }, action: ActionLeaf::Insert, hash: Some(MNT4753Fr::from_str("2").unwrap()) });
        leaves_to_process.push(OperationLeaf { coord: Coord { height: 0, idx: 16 }, action: ActionLeaf::Remove, hash: Some(MNT4753Fr::from_str("3").unwrap()) });
        leaves_to_process.push(OperationLeaf { coord: Coord { height: 0, idx: 29 }, action: ActionLeaf::Insert, hash: Some(MNT4753Fr::from_str("3").unwrap()) });
        leaves_to_process.push(OperationLeaf { coord: Coord { height: 0, idx: 16 }, action: ActionLeaf::Remove, hash: Some(MNT4753Fr::from_str("3").unwrap()) });

        let mut smt = MNT4PoseidonSMTLazy::new(
            TEST_HEIGHT_1,
            false,
            String::from("./db_leaves_mnt4"),
        ).unwrap();
        smt.process_leaves(leaves_to_process.as_slice());

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
        let mut tree = MNT4753FieldBasedMerkleTree::new(TEST_HEIGHT_1);
        tree.append(&leaves).unwrap();

        assert_eq!(tree.root(), smt.state.root, "Roots are not equal");

    }

    #[test]
    fn process_leaves_mnt6() {

        let mut leaves_to_process: Vec<OperationLeaf<MNT6753Fr>> = Vec::new();

        leaves_to_process.push(OperationLeaf { coord: Coord { height: 0, idx: 0 }, action: ActionLeaf::Insert, hash: Some(MNT6753Fr::from_str("1").unwrap()) });
        leaves_to_process.push(OperationLeaf { coord: Coord { height: 0, idx: 9 }, action: ActionLeaf::Insert, hash: Some(MNT6753Fr::from_str("2").unwrap()) });
        leaves_to_process.push(OperationLeaf { coord: Coord { height: 0, idx: 16 }, action: ActionLeaf::Remove, hash: Some(MNT6753Fr::from_str("3").unwrap()) });
        leaves_to_process.push(OperationLeaf { coord: Coord { height: 0, idx: 29 }, action: ActionLeaf::Insert, hash: Some(MNT6753Fr::from_str("3").unwrap()) });
        leaves_to_process.push(OperationLeaf { coord: Coord { height: 0, idx: 16 }, action: ActionLeaf::Remove, hash: Some(MNT6753Fr::from_str("3").unwrap()) });

        let mut smt = MNT6PoseidonSMTLazy::new(
            TEST_HEIGHT_1,
            false,
            String::from("./db_leaves_mnt6"),
        ).unwrap();
        smt.process_leaves(leaves_to_process.as_slice());

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
        let mut tree = MNT6753FieldBasedMerkleTree::new(TEST_HEIGHT_1);
        tree.append(&leaves).unwrap();

        assert_eq!(tree.root(), smt.state.root, "Roots are not equal");
    }

    #[test]
    fn test_persistency_lazy() {
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
                TEST_HEIGHT_1,
                true,
                String::from("./db_leaves_persistency_test_info_lazy"),
            ).unwrap();

            //Insert some leaves in the tree
            let mut leaves_to_process: Vec<OperationLeaf<MNT4753Fr>> = Vec::new();

            leaves_to_process.push(OperationLeaf { coord: Coord { height: 0, idx: 0 }, action: ActionLeaf::Insert, hash: Some(MNT4753Fr::from_str("1").unwrap()) });
            leaves_to_process.push(OperationLeaf { coord: Coord { height: 0, idx: 9 }, action: ActionLeaf::Insert, hash: Some(MNT4753Fr::from_str("2").unwrap()) });

            smt.process_leaves(leaves_to_process.as_slice());

            // smt gets dropped but its info should be saved
        }
        // files and directories should have been created
        assert!(Path::new("./db_leaves_persistency_test_info_lazy").exists());

        // create a non-persistent smt in another scope by restoring the previous one
        {
            let mut smt = MNT4PoseidonSMTLazy::load(
                TEST_HEIGHT_1,
                false,
                String::from("./db_leaves_persistency_test_info_lazy"),
            ).unwrap();

            // insert other leaves
            let mut leaves_to_process: Vec<OperationLeaf<MNT4753Fr>> = Vec::new();

            leaves_to_process.push(OperationLeaf { coord: Coord { height: 0, idx: 16 }, action: ActionLeaf::Insert, hash: Some(MNT4753Fr::from_str("10").unwrap()) });
            leaves_to_process.push(OperationLeaf { coord: Coord { height: 0, idx: 29 }, action: ActionLeaf::Insert, hash: Some(MNT4753Fr::from_str("3").unwrap()) });

            smt.process_leaves(leaves_to_process.as_slice());

            // if truly state has been kept, then the equality below must pass, since `root` was
            // computed in one go with another smt
            assert_eq!(root, smt.state.root);

            // smt gets dropped and info on disk destroyed
        }

        // files and directories should have been deleted
        assert!(!Path::new("./db_leaves_persistency_test_info_lazy").exists());
        assert!(!Path::new("./db_leaves_persistency_test_info_lazy_cache").exists());

        // assert being unable to restore
        assert!(MNT4PoseidonSMTLazy::load(
            TEST_HEIGHT_1,
            false,
            String::from("./db_leaves_persistency_test_info_lazy"),
        ).is_err());
    }

    #[test]
    fn test_persistency_non_lazy() {
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
            let mut smt = MNT4PoseidonSMT::new(
                TEST_HEIGHT_1,
                true,
                String::from("./db_leaves_persistency_test_info"),
            ).unwrap();

            //Insert some leaves in the tree
            smt.insert_leaf(Coord{height:0, idx:0}, MNT4753Fr::from_str("1").unwrap());
            smt.insert_leaf(Coord{height:0, idx:9}, MNT4753Fr::from_str("2").unwrap());

            // smt gets dropped but its info should be saved
        }
        // files and directories should have been created
        assert!(Path::new("./db_leaves_persistency_test_info").exists());

        // create a non-persistent smt in another scope by restoring the previous one
        {
            let mut smt = MNT4PoseidonSMT::load_batch::<MNT4753FieldBasedMerkleTreeParams>(
                TEST_HEIGHT_1,
                false,
                String::from("./db_leaves_persistency_test_info"),
            ).unwrap();

            // insert other leaves
            smt.insert_leaf(Coord { height: 0, idx: 16 }, MNT4753Fr::from_str("10").unwrap());
            smt.insert_leaf(Coord { height: 0, idx: 29 }, MNT4753Fr::from_str("3").unwrap());

            // if truly state has been kept, then the equality below must pass, since `root` was
            // computed in one go with another smt
            assert_eq!(root, smt.get_root());

            // smt gets dropped and the info on disk destroyed
        }

        // files and directories should have been deleted
        assert!(!Path::new("./db_leaves_persistency_test_info").exists());

        // assert being unable to restore
        assert!(MNT4PoseidonSMT::load_batch::<MNT4753FieldBasedMerkleTreeParams>(
            TEST_HEIGHT_1,
            false,
            String::from("./db_leaves_persistency_test_info"),
        ).is_err());
    }

    #[test]
    fn merkle_tree_path_test_mnt4() {

        let num_leaves = 32;
        let mut leaves = Vec::with_capacity(num_leaves);
        let mut leaves_for_lazy_smt = Vec::with_capacity(num_leaves);
        let mut rng = XorShiftRng::seed_from_u64(1231275789u64);

        let mut smt = MNT4PoseidonSMTLazy::new(
            TEST_HEIGHT_1,
            false,
            String::from("./db_leaves_merkle_tree_path_test_mnt4_lazy"),
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
        let mut naive_tree = MNT4753FieldBasedMerkleTree::new(TEST_HEIGHT_1);
        naive_tree.append(&leaves).unwrap();
        let root = smt.process_leaves(leaves_for_lazy_smt.as_slice());
        let naive_root = naive_tree.root();
        assert_eq!(root, naive_root);

        for i in 0..num_leaves {

            // Create and verify a FieldBasedMHTPath
            let path = smt.get_merkle_path(Coord{height: 0, idx: i});
            assert!(path.verify(smt.height(), &leaves[i], &root).unwrap());

            // Create and verify a Naive path
            let naive_path = naive_tree.generate_proof(i, &leaves[i]).unwrap();
            assert!(naive_path.verify(naive_tree.height(), &leaves[i], &naive_root ).unwrap());

            // Assert the two paths are equal
            assert_eq!(naive_path, path);

            // Check leaf index is the correct one
            assert_eq!(i, path.leaf_index());

            if i == 0 { assert!(path.is_leftmost()); } // leftmost check
            else if i == 2usize.pow(TEST_HEIGHT_1 as u32) - 1 { assert!(path.is_rightmost()) }  //rightmost check
            else { assert!(!path.is_leftmost()); assert!(!path.is_rightmost()); } // other cases check

            // Serialization/deserialization test
            let path_serialized = to_bytes!(path).unwrap();
            let path_deserialized = MNT4MerklePath::read(path_serialized.as_slice()).unwrap();
            assert_eq!(path, path_deserialized);
        }
    }

    #[test]
    fn merkle_tree_path_test_mnt6() {

        let num_leaves = 32;
        let mut leaves = Vec::with_capacity(num_leaves);
        let mut leaves_for_lazy_smt = Vec::with_capacity(num_leaves);
        let mut rng = XorShiftRng::seed_from_u64(1231275789u64);

        let mut smt = MNT6PoseidonSMTLazy::new(
            TEST_HEIGHT_1,
            false,
            String::from("./db_leaves_merkle_tree_path_test_mnt6_lazy"),
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
        let mut naive_tree = MNT6753FieldBasedMerkleTree::new(TEST_HEIGHT_1);
        naive_tree.append(&leaves).unwrap();
        let root = smt.process_leaves(leaves_for_lazy_smt.as_slice());
        let naive_root = naive_tree.root();
        assert_eq!(root, naive_root);

        for i in 0..num_leaves {

            // Create and verify a FieldBasedMHTPath
            let path = smt.get_merkle_path(Coord{height: 0, idx: i});
            assert!(path.verify(smt.height(), &leaves[i], &root).unwrap());

            // Create and verify a Naive path
            let naive_path = naive_tree.generate_proof(i, &leaves[i]).unwrap();
            assert!(naive_path.verify(naive_tree.height(), &leaves[i], &naive_root ).unwrap());

            // Assert the two paths are equal
            assert_eq!(naive_path, path);

            // Check leaf index is the correct one
            assert_eq!(i, path.leaf_index());

            if i == 0 { assert!(path.is_leftmost()); } // leftmost check
            else if i == 2usize.pow(TEST_HEIGHT_1 as u32) - 1 { assert!(path.is_rightmost()) }  //rightmost check
            else { assert!(!path.is_leftmost()); assert!(!path.is_rightmost()); } // other cases check

            // Serialization/deserialization test
            let path_serialized = to_bytes!(path).unwrap();
            let path_deserialized = MNT6MerklePath::read(path_serialized.as_slice()).unwrap();
            assert_eq!(path, path_deserialized);
        }
    }

    const TEST_HEIGHT_2: usize = 23;

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
                TEST_HEIGHT_2,
                false,
                String::from("./db_leaves_mnt4_comp_1"),
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
                TEST_HEIGHT_2,
                false,
                String::from("./db_leaves_mnt4_comp_2"),
            ).unwrap();
            let leaves_to_process2 = leaves_to_insert.clone();
            let now = Instant::now();
            root2 = smt2.process_leaves(leaves_to_process2.as_slice());
            let new_now = Instant::now();

            let duration_fast = new_now.duration_since(now).as_millis();

            println!("duration fast = {}", duration_fast);

            let leaves_to_process4 = leaves_to_remove.clone();
            let now = Instant::now();
            root4 = smt2.process_leaves(leaves_to_process4.as_slice());
            let new_now = Instant::now();

            let duration_fast = new_now.duration_since(now).as_millis();

            println!("duration fast = {}",duration_fast);

        }

        assert_eq!(root1, root2, "Roots are not equal");
        assert_eq!(root3, root4, "Roots are not equal");

        assert_eq!(root3, MNT4753_MHT_POSEIDON_PARAMETERS.nodes[23], "Sequence of roots not equal");

    }

    #[test]
    fn process_leaves_mnt6_comp() {

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
                TEST_HEIGHT_2,
                false,
                String::from("./db_leaves_mnt6_comp_1"),
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
                TEST_HEIGHT_2,
                false,
                String::from("./db_leaves_mnt6_comp_2"),
            ).unwrap();
            let leaves_to_process2 = leaves_to_insert.clone();
            let now = Instant::now();
            root2 = smt2.process_leaves(leaves_to_process2.as_slice());
            let new_now = Instant::now();

            let duration_fast = new_now.duration_since(now).as_millis();

            println!("duration fast = {}", duration_fast);

            let leaves_to_process4 = leaves_to_remove.clone();
            let now = Instant::now();
            root4 = smt2.process_leaves(leaves_to_process4.as_slice());
            let new_now = Instant::now();

            let duration_fast = new_now.duration_since(now).as_millis();

            println!("duration fast = {}",duration_fast);

        }

        assert_eq!(root1, root2, "Roots are not equal");
        assert_eq!(root3, root4, "Roots are not equal");

        assert_eq!(root3, MNT6753_MHT_POSEIDON_PARAMETERS.nodes[23], "Sequence of roots not equal");

    }
}