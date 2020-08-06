use crate::merkle_tree::field_based_mht::smt::{SmtPoseidonParameters, Coord, OperationLeaf, BigMerkleTreeState};
use crate::{PoseidonParameters, BatchFieldBasedHash, merkle_tree};
use crate::crh::poseidon::batched_crh::PoseidonBatchHash;
use crate::merkle_tree::field_based_mht::smt::ActionLeaf::Remove;
use crate::merkle_tree::field_based_mht::smt::parameters::{MNT4753SmtPoseidonParameters, MNT6753SmtPoseidonParameters};
use crate::crh::poseidon::parameters::{MNT4753PoseidonParameters, MNT6753PoseidonParameters};

use rocksdb::{DB, Options};

use std::collections::HashSet;
use std::marker::PhantomData;

use std::fs;

use algebra::{PrimeField, MulShort};
use algebra::{ToBytes, to_bytes, FromBytes};
use algebra::fields::mnt6753::Fr as MNT6753Fr;
use algebra::fields::mnt4753::Fr as MNT4753Fr;
use crate::merkle_tree::field_based_mht::smt::big_merkle_tree::BigMerkleTree;
use crate::merkle_tree::field_based_mht::smt::error::Error;

#[derive(Debug)]
pub struct LazyBigMerkleTree<F: PrimeField, T: SmtPoseidonParameters<Fr=F>, P: PoseidonParameters<Fr=F>> {
    // if unset, all DBs and tree internal state will be deleted when an instance of this struct
    // gets dropped
    persistent: bool,
    // path to state required to restore the tree
    state_path: Option<String>,
    // tree in-memory state
    state: BigMerkleTreeState<F, T>,
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

    _field: PhantomData<F>,
    _parameters: PhantomData<T>,
    _poseidon_parameters: PhantomData<P>,
}

impl<F: PrimeField, T: SmtPoseidonParameters<Fr=F>, P: PoseidonParameters<Fr=F>> Drop for LazyBigMerkleTree<F, T, P> {
    fn drop(&mut self) {

        if !self.persistent {

            if self.state_path.is_some() {
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
}

impl<F: PrimeField + MulShort, T: SmtPoseidonParameters<Fr=F>, P: PoseidonParameters<Fr=F>> LazyBigMerkleTree<F, T, P> {
    // Creates a new tree of specified `width`.
    // If `persistent` is specified, then DBs will be kept on disk and the tree state will be saved
    // so that the tree can be restored any moment later. Otherwise, no state will be saved on file
    // and the DBs will be deleted.
    pub fn new_unitialized(width: usize, persistent: bool, state_path: Option<String>, path_db: String, path_cache: String) -> Result<Self, Error> {

        // If the tree must be persistent, that a path to which save the tree state must be
        // specified.
        if !persistent { assert!(state_path.is_none()) } else { assert!(state_path.is_some()) }

        let height = width as f64;
        let height = height.log(T::MERKLE_ARITY as f64) as usize;
        let state = BigMerkleTreeState::<F, T>::get_default_state(width, height);
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
            _field: PhantomData,
            _parameters: PhantomData,
            _poseidon_parameters: PhantomData,
        })
    }

    // Creates a new tree starting from state at `state_path` and DBs at `path_db` and `path_cache`.
    // The new tree may be persistent or not, the actions taken in both cases are the same as
    // in `new_unitialized()`.
    pub fn new(persistent: bool, state_path: String, path_db: String, path_cache: String) -> Result<Self, Error> {

        let state = {
            let state_file = fs::File::open(state_path.clone())
                .map_err(|e| Error::Other(e.to_string()))?;
            BigMerkleTreeState::<F, T>::read(state_file)
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
            _field: PhantomData,
            _parameters: PhantomData,
            _poseidon_parameters: PhantomData,
        })
    }

    pub fn height(&self) -> usize { self.height }

    pub fn insert_to_cache(&self, coord: Coord, data:F) {
        let elem = to_bytes!(data).unwrap();
        let index = bincode::serialize(&coord).unwrap();
        self.database.put(index, elem).unwrap();
    }

    pub fn contains_key_in_cache(&self, coord:Coord) -> bool {
        let coordinates = bincode::serialize(&coord).unwrap();
        match self.database.get(coordinates) {
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

    pub fn get_from_cache(&self, coord:Coord) -> Option<F> {
        let coordinates = bincode::serialize(&coord).unwrap();
        match self.database.get(coordinates) {
            Ok(Some(value)) => {
                let retrieved_elem = F::read(value.as_slice()).unwrap();
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

    pub fn remove_from_cache(&self, coord: Coord) -> Option<F>{
        let coordinates = bincode::serialize(&coord).unwrap();
        match self.database.get(coordinates.clone()) {
            Ok(Some(value)) => {
                let retrieved_elem = F::read(value.as_slice()).unwrap();
                let res = self.database.delete(coordinates.clone());
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

    pub fn insert_to_db(&self, idx: usize, data: F) {
        let elem = to_bytes!(data).unwrap();
        let index = bincode::serialize(&idx).unwrap();
        self.database.put(index, elem).unwrap();
    }

    pub fn get_from_db(&self, idx: usize) -> Option<F>{
        let index = bincode::serialize(&idx).unwrap();
        match self.database.get(index) {
            Ok(Some(value)) => {
                let retrieved_elem = F::read(value.as_slice()).unwrap();
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

    pub fn remove_from_db(&self, idx: usize) -> Option<F>{
        let index = bincode::serialize(&idx).unwrap();
        match self.database.get(index.clone()) {
            Ok(Some(value)) => {
                let retrieved_elem = F::read(value.as_slice()).unwrap();
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
    pub fn node(&mut self, coord: Coord) -> F {

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
                let left_hash: F;
                if let Some(i) = left_child {
                    left_hash = i;
                } else {
                    left_hash = T::EMPTY_HASH_CST[0];
                }

                let right_child_idx = left_child_idx + 1;
                let right_child = self.get_from_db(right_child_idx);
                let right_hash: F;
                if let Some(i) = right_child {
                    right_hash = i;
                } else {
                    right_hash = T::EMPTY_HASH_CST[0];
                }
                node_hash = merkle_tree::field_based_mht::smt::big_merkle_tree::BigMerkleTree::<F, T, P>::poseidon_hash(left_hash, right_hash);
            } else {
                let height_child = coord.height - 1;
                let left_child_idx = coord.idx * T::MERKLE_ARITY;
                let coord_left = Coord { height: height_child, idx: left_child_idx };
                let left_child_hash = LazyBigMerkleTree::node(self, coord_left);

                let right_child_idx = left_child_idx + 1;
                let coord_right = Coord { height: height_child, idx: right_child_idx };
                let right_child_hash = LazyBigMerkleTree::node(self, coord_right);

                node_hash = merkle_tree::field_based_mht::smt::big_merkle_tree::BigMerkleTree::<F, T, P>::poseidon_hash(left_child_hash, right_child_hash);
            }
            return node_hash;
        }

        res.unwrap()
    }

    pub fn get_root(&self) -> F {
        self.state.root.clone()
    }

    pub fn remove_node_from_cache(&mut self, coord: Coord) {
        self.remove_from_cache(coord);
    }

    pub fn batch_poseidon_hash(input: Vec<F>) -> Vec<F> {
        let output_vec = PoseidonBatchHash::<F, P>::batch_evaluate(&input);
        output_vec.unwrap()
    }

    pub fn process_leaves (&mut self, vec_leaf_op: Vec<OperationLeaf<F>>) -> F {

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
                both_children_present.push(true);
            } else {
                both_children_present.push(false);
            }
        }
        // Process the input_vec using batch Poseidon hash
        let output_vec = merkle_tree::field_based_mht::smt::big_lazy_merkle_tree::LazyBigMerkleTree::<F, T, P>::batch_poseidon_hash(input_vec);
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
                let left_hash: F;
                let right_hash: F;
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
            let output_vec = merkle_tree::field_based_mht::smt::big_lazy_merkle_tree::LazyBigMerkleTree::<F, T, P>::batch_poseidon_hash(input_vec);

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
        self.state.root
    }

}
pub type MNT4PoseidonSmt = BigMerkleTree<MNT4753Fr, MNT4753SmtPoseidonParameters, MNT4753PoseidonParameters>;
pub type MNT4PoseidonSmtLazy = LazyBigMerkleTree<MNT4753Fr, MNT4753SmtPoseidonParameters, MNT4753PoseidonParameters>;
pub type MNT6PoseidonSmt = BigMerkleTree<MNT6753Fr, MNT6753SmtPoseidonParameters, MNT6753PoseidonParameters>;
pub type MNT6PoseidonSmtLazy = LazyBigMerkleTree<MNT6753Fr, MNT6753SmtPoseidonParameters, MNT6753PoseidonParameters>;

#[cfg(test)]
mod test {
    use crate::merkle_tree::field_based_mht::smt::{MNT4PoseidonHash, OperationLeaf, Coord, ActionLeaf, SmtPoseidonParameters, MNT6PoseidonHash};
    use crate::merkle_tree::field_based_mht::{FieldBasedMerkleTreeConfig, FieldBasedMerkleHashTree};
    use crate::merkle_tree::field_based_mht::smt::parameters::{MNT4753SmtPoseidonParameters, MNT6753SmtPoseidonParameters};
    use crate::merkle_tree::field_based_mht::smt::big_lazy_merkle_tree::{MNT4PoseidonSmt, MNT4PoseidonSmtLazy, MNT6PoseidonSmt, MNT6PoseidonSmtLazy};

    use algebra::fields::mnt6753::Fr as MNT6753Fr;
    use algebra::fields::mnt4753::Fr as MNT4753Fr;
    use std::str::FromStr;
    use std::time::Instant;

    use rand_xorshift::XorShiftRng;
    use rand::{SeedableRng, RngCore};
    use rand::rngs::OsRng;

    struct MNT4753FieldBasedMerkleTreeParams;
    struct MNT6753FieldBasedMerkleTreeParams;

    impl FieldBasedMerkleTreeConfig for MNT4753FieldBasedMerkleTreeParams {
        const HEIGHT: usize = 6;
        type H = MNT4PoseidonHash;
    }

    impl FieldBasedMerkleTreeConfig for MNT6753FieldBasedMerkleTreeParams {
        const HEIGHT: usize = 6;
        type H = MNT6PoseidonHash;
    }

    type MNT4753FieldBasedMerkleTree = FieldBasedMerkleHashTree<MNT4753FieldBasedMerkleTreeParams>;
    type MNT6753FieldBasedMerkleTree = FieldBasedMerkleHashTree<MNT6753FieldBasedMerkleTreeParams>;

    #[test]
    fn process_leaves_mnt4_comp() {
        use algebra::{
            fields::mnt4753::Fr,
            UniformRand,
        };

        let num_leaves = 2usize.pow(23);
        let mut rng1 = XorShiftRng::seed_from_u64(9174123u64);
        let mut leaves_to_insert: Vec<OperationLeaf<MNT4753Fr>> = Vec::new();
        let mut leaves_to_remove: Vec<OperationLeaf<MNT4753Fr>> = Vec::new();

        let n = 1000;

        for _i in 0..n {
            let random: u64 = OsRng.next_u64();
            let idx = random % num_leaves as u64;
            let elem = Fr::rand(&mut rng1);

            leaves_to_insert.push(OperationLeaf { coord: Coord { height: 0, idx: idx as usize }, action: ActionLeaf::Insert, hash: Some(elem.clone()) });
            leaves_to_remove.push(OperationLeaf { coord: Coord { height: 0, idx: idx as usize }, action: ActionLeaf::Remove, hash: None });
        }

        // Insertion

        let root1;
        let root2;
        let root3;
        let root4;
        {
            let mut smt1 = MNT4PoseidonSmt::new_unitialized(
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
            let mut smt2 = MNT4PoseidonSmtLazy::new_unitialized(
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

        assert_eq!(root3, MNT4753SmtPoseidonParameters::EMPTY_HASH_CST[23], "Sequence of roots not equal");

    }

    #[test]
    fn process_leaves_mnt4() {
        use algebra::{
            fields::mnt4753::Fr, Field,
        };

        let num_leaves = 32;
        let mut leaves_to_process: Vec<OperationLeaf<MNT4753Fr>> = Vec::new();

        leaves_to_process.push(OperationLeaf { coord: Coord { height: 0, idx: 0 }, action: ActionLeaf::Insert, hash: Some(MNT4753Fr::from_str("1").unwrap()) });
        leaves_to_process.push(OperationLeaf { coord: Coord { height: 0, idx: 9 }, action: ActionLeaf::Insert, hash: Some(MNT4753Fr::from_str("2").unwrap()) });
        leaves_to_process.push(OperationLeaf { coord: Coord { height: 0, idx: 16 }, action: ActionLeaf::Remove, hash: Some(MNT4753Fr::from_str("3").unwrap()) });
        leaves_to_process.push(OperationLeaf { coord: Coord { height: 0, idx: 29 }, action: ActionLeaf::Insert, hash: Some(MNT4753Fr::from_str("3").unwrap()) });
        leaves_to_process.push(OperationLeaf { coord: Coord { height: 0, idx: 16 }, action: ActionLeaf::Remove, hash: Some(MNT4753Fr::from_str("3").unwrap()) });

        let mut smt = MNT4PoseidonSmtLazy::new_unitialized(
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
            let f = Fr::zero();
            leaves.push(f);
        }
        leaves.push(MNT4753Fr::from_str("2").unwrap());
        for _ in 10..29 {
            let f = Fr::zero();
            leaves.push(f);
        }
        leaves.push(MNT4753Fr::from_str("3").unwrap());
        for _ in 30..31 {
            let f = Fr::zero();
            leaves.push(f);
        }
        let tree = MNT4753FieldBasedMerkleTree::new(&leaves).unwrap();

        assert_eq!(tree.root.unwrap(), smt.state.root, "Roots are not equal");

    }

    #[test]
    fn process_leaves_mnt6() {
        use algebra::{
            fields::mnt6753::Fr,
            UniformRand,
        };

        let num_leaves = 2usize.pow(23);
        let mut rng1 = XorShiftRng::seed_from_u64(9174123u64);
        let mut leaves_to_insert: Vec<OperationLeaf<MNT6753Fr>> = Vec::new();
        let mut leaves_to_remove: Vec<OperationLeaf<MNT6753Fr>> = Vec::new();

        let n = 1000;

        for _i in 0..n {
            let random: u64 = OsRng.next_u64();
            let idx = random % num_leaves as u64;
            let elem = Fr::rand(&mut rng1);

            leaves_to_insert.push(OperationLeaf { coord: Coord { height: 0, idx: idx as usize }, action: ActionLeaf::Insert, hash: Some(elem.clone()) });
            leaves_to_remove.push(OperationLeaf { coord: Coord { height: 0, idx: idx as usize }, action: ActionLeaf::Remove, hash: None });
        }

        // Insertion

        let root1;
        let root2;
        let root3;
        let root4;
        {
            let mut smt1 = MNT6PoseidonSmt::new_unitialized(
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
            let mut smt2 = MNT6PoseidonSmtLazy::new_unitialized(
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

        assert_eq!(root3, MNT6753SmtPoseidonParameters::EMPTY_HASH_CST[23], "Sequence of roots not equal");

    }

    #[test]
    fn process_leaves_mnt6_comp() {
        use algebra::{
            fields::mnt6753::Fr, Field,
        };

        let num_leaves = 32;
        let mut leaves_to_process: Vec<OperationLeaf<MNT6753Fr>> = Vec::new();

        leaves_to_process.push(OperationLeaf { coord: Coord { height: 0, idx: 0 }, action: ActionLeaf::Insert, hash: Some(MNT6753Fr::from_str("1").unwrap()) });
        leaves_to_process.push(OperationLeaf { coord: Coord { height: 0, idx: 9 }, action: ActionLeaf::Insert, hash: Some(MNT6753Fr::from_str("2").unwrap()) });
        leaves_to_process.push(OperationLeaf { coord: Coord { height: 0, idx: 16 }, action: ActionLeaf::Remove, hash: Some(MNT6753Fr::from_str("3").unwrap()) });
        leaves_to_process.push(OperationLeaf { coord: Coord { height: 0, idx: 29 }, action: ActionLeaf::Insert, hash: Some(MNT6753Fr::from_str("3").unwrap()) });
        leaves_to_process.push(OperationLeaf { coord: Coord { height: 0, idx: 16 }, action: ActionLeaf::Remove, hash: Some(MNT6753Fr::from_str("3").unwrap()) });

        let mut smt = MNT6PoseidonSmtLazy::new_unitialized(
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
            let f = Fr::zero();
            leaves.push(f);
        }
        leaves.push(MNT6753Fr::from_str("2").unwrap());
        for _ in 10..29 {
            let f = Fr::zero();
            leaves.push(f);
        }
        leaves.push(MNT6753Fr::from_str("3").unwrap());
        for _ in 30..31 {
            let f = Fr::zero();
            leaves.push(f);
        }
        let tree = MNT6753FieldBasedMerkleTree::new(&leaves).unwrap();

        assert_eq!(tree.root.unwrap(), smt.state.root, "Roots are not equal");
    }


}