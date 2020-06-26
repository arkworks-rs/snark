mod parameters;

use crate::merkle_tree::field_based_mht::smt::error::Error;
use std::collections::HashMap;
use std::collections::HashSet;
use crate::{FieldBasedHashParameters, PoseidonHash, PoseidonParameters, FieldBasedHash, merkle_tree, BatchFieldBasedHash};

use rocksdb::DB;

use algebra::fields::mnt6753::Fr as MNT6753Fr;
use algebra::fields::mnt4753::Fr as MNT4753Fr;

use algebra::{ToBytes, to_bytes};

use algebra::{PrimeField, MulShort};
use std::marker::PhantomData;
use crate::crh::poseidon::parameters::{MNT4753PoseidonParameters, MNT6753PoseidonParameters};
use crate::crh::poseidon::batched_crh::PoseidonBatchHash;
use crate::merkle_tree::field_based_mht::smt::ActionLeaf::{Remove, Insert};

use crate::merkle_tree::field_based_mht::smt::parameters::{MNT4753SmtPoseidonParameters, MNT6753SmtPoseidonParameters};
use std::fs;

use serde::{Serialize,Deserialize};

pub mod error;

pub type MNT4PoseidonHash = PoseidonHash<MNT4753Fr, MNT4753PoseidonParameters>;
pub type MNT6PoseidonHash = PoseidonHash<MNT6753Fr, MNT6753PoseidonParameters>;

pub trait SmtPoseidonParameters: 'static + FieldBasedHashParameters {
    // The arity of the Sparse Merkle Tree
    const MERKLE_ARITY: usize;
    // The pre-computed hashes of the empty nodes for the different levels of the SMT
    const EMPTY_HASH_CST: &'static [Self::Fr];
}

#[derive(Debug, PartialEq, Eq, Hash, Copy, Clone)]
pub enum ActionLeaf {
    Insert,
    Remove,
}

#[derive(Debug)]
pub struct BigMerkleTree<F: PrimeField, T: SmtPoseidonParameters<Fr=F>, P: PoseidonParameters<Fr=F>> {

    // the number of leaves
    width: usize,
    // the height of the Merkle tree
    height: usize,
    // path to the db
    path_db: &'static str,
    // stores the leaves
    database: DB,
    // path to the cache
    path_cache: &'static str,
    // stores the cached nodes
    db_cache: DB,
    cache: HashMap<Coord, F>,
    // stores the nodes of the path
    cache_path: HashMap<Coord, F>,
    // indicates which nodes are present the Merkle tree
    present_node: HashSet<Coord>,
    // cache of nodes processed in parallel
    cache_parallel: HashMap<Coord, F>,
    // root of the Merkle tree
    root: F,

    // for debug purposes, print intermediate results
    print_verbose: bool,

    _field: PhantomData<F>,
    _parameters: PhantomData<T>,
    _poseidon_parameters: PhantomData<P>,
}

#[derive(Debug, PartialEq, Eq, Hash, Copy, Clone, Serialize, Deserialize)]
// Coordinates system for identifying a node
pub struct Coord {
    // height in the Merkle tree (0 -> leaves)
    height: usize,
    // the index of the node in that level
    idx: usize,
}

impl Coord {
    pub fn new(height: usize, idx: usize) -> Self {
        Self { height, idx }
    }
}

#[derive(Debug, PartialEq, Eq, Hash, Copy, Clone)]
// Action associated to the leaf
pub struct OperationLeaf <F: PrimeField>{
    coord: Coord,
    action: ActionLeaf,
    hash: Option<F>,
}

impl<F: PrimeField> OperationLeaf<F> {
    pub fn new(height: usize, idx: usize, action: ActionLeaf, hash: Option<F>) -> Self {
        Self { coord: Coord { height, idx }, action, hash}
    }
}

impl<F: PrimeField, T: SmtPoseidonParameters<Fr=F>, P: PoseidonParameters<Fr=F>> Drop for BigMerkleTree<F, T, P> {
    fn drop(&mut self) {
        // Deletes the folder containing the db
        match fs::remove_dir_all(self.path_db) {
            Ok(_) => (),
            Err(e) => {
                println!("Error deleting the folder containing the db. {}", e);
            }
        };
        // Deletes the folder containing the cache
        match fs::remove_dir_all(self.path_cache) {
            Ok(_) => (),
            Err(e) => {
                println!("Error deleting the folder containing the db. {}", e);
            }
        };
    }
}

// Assumption: MERKLE_ARITY == 2
impl<F: PrimeField + MulShort, T: SmtPoseidonParameters<Fr=F>, P: PoseidonParameters<Fr=F>> BigMerkleTree<F, T, P> {
    pub fn new(width: usize) -> Result<Self, Error> {
        let height = width as f64;
        let height = height.log(T::MERKLE_ARITY as f64) as usize;
        let path_db = "/home/mkaihara/Documents/dev/ginger-lib/primitives/src/merkle_tree/field_based_mht/smt/rocksdb_storage";
        let database = DB::open_default(path_db).unwrap();
        let path_cache = "/home/mkaihara/Documents/dev/ginger-lib/primitives/src/merkle_tree/field_based_mht/smt/rocksdb_cache_storage";
        let db_cache = DB::open_default(path_cache).unwrap();
        let cache = HashMap::new();
        let cache_path = HashMap::new();
        let present_node = HashSet::new();
        let cache_parallel = HashMap::new();
        let root = T::EMPTY_HASH_CST[height];
        // set to true to print intermediate results for debug
        let print_verbose = false;
        Ok(BigMerkleTree {
            width,
            height,
            path_db,
            database,
            path_cache,
            db_cache,
            cache,
            cache_path,
            present_node,
            cache_parallel,
            root,
            print_verbose,
            _field: PhantomData,
            _parameters: PhantomData,
            _poseidon_parameters: PhantomData,
        })
    }

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


    // Inserts a leaf in the Merkle tree. Only used for single operation
    // Updates the Merkle tree on the path from the leaf to the root
    pub fn insert_leaf(&mut self, coord: Coord, leaf: F) {

        // check that the index of the leaf to be inserted is less than the width of the Merkle tree
        assert!(coord.idx < self.width, "Leaf index out of bound.");
        // check that the coordinates of the node corresponds to the leaf level
        assert_eq!(coord.height, 0, "Coord of the node does not correspond to leaf level");

        if self.present_node.contains(&coord) {
            let old_leaf = self.get_from_db(coord.idx);
            let old_hash;
            if let Some(i) = old_leaf {
                old_hash = i;
            } else {
                old_hash = T::EMPTY_HASH_CST[0];
            }
            if old_hash != leaf {
                self.insert_to_db(coord.idx, leaf);
                if self.print_verbose { println!("insert_leaf: Leaf already present. Update leaf in db {:?}", coord) };
                self.cache_path.clear();
                self.cache_path.insert(coord, leaf);
                if self.print_verbose { println!("insert_leaf: Insert leaf in cache_path {:?}", coord); }
                BigMerkleTree::update_tree(self, coord);
            } else {
                if self.print_verbose { println!("insert_leaf: Leaf already present. Leaf equal. Does not update leaf in db {:?}", coord) };
            }
        } else {
            // mark as non empty leaf
            self.present_node.insert(coord);
            if self.print_verbose { println!("insert_leaf: New leaf. Insert leaf in non_empty set {:?}", coord); }
            self.insert_to_db(coord.idx, leaf);
            if self.print_verbose { println!("insert_leaf: Insert new leaf in db {:?}", coord); }
            self.cache_path.clear();
            self.cache_path.insert(coord, leaf);
            if self.print_verbose { println!("insert_leaf: Insert new leaf in cache_path {:?}", coord); }
            BigMerkleTree::update_tree(self, coord);
        }
    }

    // Removes a leaf in the Merkle tree. Only used for single operation
    // Updates the Merkle tree on the path from the leaf to the root
    pub fn remove_leaf(&mut self, coord: Coord) {

        // check that the index of the leaf to be inserted is less than the width of the Merkle tree
        assert!(coord.idx < self.width, "Leaf index out of bound.");
        // check that the coordinates of the node corresponds to the leaf level
        assert_eq!(coord.height, 0, "Coord of the node does not correspond to leaf level");

        // take that leaf from the non-empty set
        self.present_node.remove(&coord);
        if self.print_verbose { println!("remove_leaf: Remove leaf from non_empty set {:?}", coord); }
        self.cache_path.clear();
        self.cache_path.insert(coord, T::EMPTY_HASH_CST[0]);
        if self.print_verbose { println!("remove_leaf: Insert empty leaf in the cache_path {:?}", coord); }
        // removes the leaf from the db
        let res = self.remove_from_db(coord.idx);
        if self.print_verbose { println!("remove_leaf: Remove leaf from db {:?}", coord); }
        // if it was in the db, update the tree
        if res != None {
            BigMerkleTree::update_tree(self, coord);
        }
    }

    // Updates the tree visiting the parent nodes from the leaf to the root
    // Calculates the hash and caches it
    pub fn update_tree(&mut self, coord: Coord) {

        // Process the node of level 1 with the inserted/removed leaf
        // check whether the hash corresponds to the left or right child
        let mut idx = coord.idx;
        let left_child_coord: Coord;
        let right_child_coord: Coord;
        let left_child_idx: usize;
        let right_child_idx: usize;
        let left_child: Option<F>;
        let right_child: Option<F>;
        let left_hash: F;
        let right_hash: F;
        let mut height = 0;

        assert_eq!(T::MERKLE_ARITY,2, "Arity of the Merkle tree is not 2.");

        if self.print_verbose { println!("update_tree: ===============================================================>"); }

        if idx % T::MERKLE_ARITY == 0 {
            left_child_idx = idx;
            left_child_coord = Coord { height, idx: left_child_idx };
            // get the left child from the cache_path
            let hash = self.cache_path.get(&left_child_coord);
            if let Some(i) = hash {
                left_hash = *i;
                if self.print_verbose { println!("update_tree: Get leaf from cache_path {:?}", left_child_coord); }
            } else {
                left_hash = T::EMPTY_HASH_CST[0];
                if self.print_verbose { println!("update_tree: Set empty leaf {:?}", left_child_coord); }
            }
            right_child_idx = idx + 1;
            right_child_coord = Coord { height, idx: right_child_idx };
            right_child = self.get_from_db(right_child_idx);
            if let Some(i) = right_child {
                right_hash = i;
                if self.print_verbose { println!("update_tree: Get leaf from db {:?}", right_child_coord); }
            } else {
                right_hash = T::EMPTY_HASH_CST[0];
                if self.print_verbose { println!("update_tree: Set empty leaf {:?}", right_child_coord); }
            }
        } else {
            right_child_idx = idx;
            right_child_coord = Coord { height, idx: right_child_idx };
            // get the right child from the cache path
            let hash = self.cache_path.get(&right_child_coord);
            if let Some(i) = hash {
                right_hash = *i;
                if self.print_verbose { println!("update_tree: Get leaf from cache_path {:?}", right_child_coord); }
            } else {
                right_hash = T::EMPTY_HASH_CST[0];
                if self.print_verbose { println!("update_tree: Set empty leaf {:?}", right_child_coord); }
            }
            left_child_idx = idx - 1;
            left_child_coord = Coord { height, idx: left_child_idx };
            left_child = self.get_from_db(left_child_idx);
            if let Some(i) = left_child {
                left_hash = i;
                if self.print_verbose { println!("update_tree: Get leaf from cache_path {:?}", left_child_coord); }
            } else {
                left_hash = T::EMPTY_HASH_CST[0];
                if self.print_verbose { println!("update_tree: Set empty leaf {:?}", left_child_coord); }
            }
        }

        height += 1;
        idx = idx / T::MERKLE_ARITY;
        let parent_coord = Coord { height, idx };

        let mut node_hash: F;
        if (!self.present_node.contains(&left_child_coord)) & (!self.present_node.contains(&right_child_coord)) {
            node_hash = T::EMPTY_HASH_CST[height];
            if self.print_verbose { println!("update_tree: -----> Visit parent node {:?}. Node is empty.", parent_coord); }

            // insert the parent node into the cache_path
            self.cache_path.insert(parent_coord, node_hash);
            if self.print_verbose { println!("update_tree: insert empty node in cache_path {:?}", parent_coord); }

            // Both children are empty leaves
            // remove the parent node from the presence set
            self.present_node.remove(&parent_coord.clone());
            if self.print_verbose { println!("update_tree: remove node from empty_list {:?}", parent_coord); }
            // remove node from cache
            BigMerkleTree::remove_node_from_cache(self, parent_coord.clone());
            if self.print_verbose { println!("update_tree: remove node from cache {:?}", parent_coord); }
        } else {

            if self.cache_parallel.contains_key(&parent_coord) {
                node_hash = *self.cache_parallel.get(&parent_coord).unwrap();
                self.cache_parallel.remove(&parent_coord);
            } else {
                // compute the hash of the node with the hashes of the children
                node_hash = merkle_tree::field_based_mht::smt::BigMerkleTree::<F, T, P>::poseidon_hash(left_hash, right_hash);
            }
            if self.print_verbose { println!("update_tree: -----> Visit parent node {:?}. Node is not empty.", parent_coord); }
            // insert the parent node into the cache_path
            self.cache_path.insert(parent_coord, node_hash);
            if self.print_verbose { println!("insert_leaf: insert node in cache_path {:?}", parent_coord); }
            // set the parent as present
            self.present_node.insert(parent_coord.clone());
            if self.print_verbose { println!("update_tree: insert node in presence set {:?}", parent_coord); }
            // Both children are not empty leaves
            if (self.present_node.contains(&left_child_coord)) & (self.present_node.contains(&right_child_coord)) {
                // cache the node
                BigMerkleTree::insert_node_in_cache(self, parent_coord.clone(), node_hash);
                if self.print_verbose { println!("update_tree: insert node in cache {:?}", parent_coord); }
            }
        }

        // Process level >= 2
        while height != self.height {
            // go to parent node
            height += 1;
            let idx_child = idx;
            idx = idx / T::MERKLE_ARITY;
            let parent_coord = Coord { height, idx };
            if self.print_verbose { println!("update_tree: -----> Visit parent node {:?}. ", parent_coord); }

            // retrieve the left child
            let left_child_idx = parent_coord.idx * T::MERKLE_ARITY;
            let left_child_height = parent_coord.height - 1;
            let left_child_coord = Coord { height: left_child_height, idx: left_child_idx };
            let left_hash: F;
            if left_child_idx == idx_child {
                // It was cached in the previous iteration. Get it from the cache_path
                left_hash = *self.cache_path.get(&left_child_coord).unwrap();
                if self.print_verbose { println!("update_tree: Get leaf from cache_path {:?}. ", left_child_coord); }
            } else {
                if self.print_verbose { println!("update_tree: Get leaf calling node {:?}. ", left_child_coord); }
                left_hash = BigMerkleTree::node(self, left_child_coord);
            }

            // retrieve the right child
            let right_child_idx = left_child_idx + 1;
            let right_child_height = left_child_height;
            let right_child_coord = Coord { height: right_child_height, idx: right_child_idx };
            let right_hash: F;
            if right_child_idx == idx_child {
                // It was cached in the previous iteration. Get it from the cache_path
                right_hash = *self.cache_path.get(&right_child_coord).unwrap();
                if self.print_verbose { println!("update_tree: Get leaf from cache_path {:?}. ", right_child_coord); }
            } else {
                if self.print_verbose { println!("update_tree: Get leaf calling node {:?}. ", right_child_coord); }
                right_hash = BigMerkleTree::node(self, right_child_coord);
            }

            if (!self.present_node.contains(&left_child_coord)) & (!self.present_node.contains(&right_child_coord)) {
                node_hash = T::EMPTY_HASH_CST[height];
                if self.print_verbose { println!("update_tree: Node is empty {:?}. ", parent_coord); }
                // insert the parent node into the cache_path
                self.cache_path.insert(parent_coord, node_hash);
                if self.print_verbose { println!("update_tree: Insert empty node in cache_path {:?}. ", parent_coord); }
                // remove node from non_empty set
                self.present_node.remove(&parent_coord.clone());
                if self.print_verbose { println!("update_tree: Remove node from presence set {:?}. ", parent_coord); }
            } else {
                // compute the hash of the parent node based on the hashes of the children
                node_hash = merkle_tree::field_based_mht::smt::BigMerkleTree::<F, T, P>::poseidon_hash(left_hash, right_hash);
                if self.print_verbose { println!("update_tree: Hash of the Node is computed {:?}. ", parent_coord); }
                // insert the parent node into the cache_path
                self.cache_path.insert(parent_coord, node_hash);
                if self.print_verbose { println!("update_tree: Insert node in the cache_path {:?}. ", parent_coord); }
                // both children are present
                if self.present_node.contains(&left_child_coord) & self.present_node.contains(&right_child_coord) {
                    // set the parent node as non_empty
                    self.present_node.insert(parent_coord.clone());
                    if self.print_verbose { println!("update_tree: Insert parent node in presence set {:?}", parent_coord); }
                    // children not empty leaves, then cache the parent node
                    BigMerkleTree::insert_node_in_cache(self, parent_coord.clone(), node_hash);
                    if self.print_verbose { println!("update_tree: Insert parent node in cache {:?}", parent_coord); }
                    // cache the children
                    BigMerkleTree::insert_node_in_cache(self, left_child_coord.clone(), left_hash);
                    if self.print_verbose { println!("update_tree: Insert left child node in cache {:?}", left_child_coord); }
                    BigMerkleTree::insert_node_in_cache(self, right_child_coord.clone(), right_hash);
                    if self.print_verbose { println!("update_tree: Insert right child node in cache {:?}", right_child_coord); }
                } else {
                    // either child not equal to empty leaves, include the parent node in a non-empty set
                    self.present_node.insert(parent_coord.clone());
                    if self.print_verbose { println!("update_tree: Insert parent node in presence set {:?}", parent_coord); }
                }
            }
            if (!self.present_node.contains(&left_child_coord)) | (!self.present_node.contains(&right_child_coord)) {
                if self.contains_key_in_cache(parent_coord) {
                    // remove subtree from cache
                    if self.print_verbose { println!("update_tree: Remove superfluous node from subtree starting at node {:?}", parent_coord); }
                    BigMerkleTree::remove_subtree_from_cache(self, parent_coord);
                }
            }
        }
        self.root = node_hash;
    }

    pub fn remove_subtree_from_cache(&mut self, coord: Coord) {

        assert_eq!(T::MERKLE_ARITY,2, "Arity of the Merkle tree is not 2.");

        if coord.height == 1 {
            let left_child_idx = coord.idx * T::MERKLE_ARITY;
            let left_child_height = 0;
            let left_coord = Coord { height: left_child_height, idx: left_child_idx };

            let right_child_idx = left_child_idx + 1;
            let right_child_height = 0;
            let right_coord = Coord { height: right_child_height, idx: right_child_idx };

            if (!self.present_node.contains(&left_coord)) | (!self.present_node.contains(&right_coord)) {
                if self.contains_key_in_cache(coord) {
                    BigMerkleTree::remove_node_from_cache(self, coord);
                    if self.print_verbose { println!("remove_subtree_from_cache: Remove node from cache {:?}", coord); }
                }
            }
            return;
        }
        if coord.height > 1 {
            let left_child_idx = coord.idx * T::MERKLE_ARITY;
            let left_child_height = coord.height - 1;
            let left_coord = Coord { height: left_child_height, idx: left_child_idx };

            let right_child_idx = left_child_idx + 1;
            let right_child_height = left_child_height;
            let right_coord = Coord { height: right_child_height, idx: right_child_idx };

            if (!self.present_node.contains(&left_coord)) | (!self.present_node.contains(&right_coord)) {
                if self.contains_key_in_cache(coord) {
                    BigMerkleTree::remove_node_from_cache(self, coord);
                    if self.print_verbose { println!("remove_subtree_from_cache: Remove node from cache {:?}", coord); }

                    if self.present_node.contains(&left_coord) {
                        if left_coord.height > 0 {

                        }
                        BigMerkleTree::remove_subtree_from_cache(self, left_coord);
                    }
                    if self.present_node.contains(&right_coord) {
                        BigMerkleTree::remove_subtree_from_cache(self, right_coord);
                    }

                }
            }

        }
    }

    pub fn check_b_plus_caching_level_down(&mut self, coord: Coord) {

        let left_child_idx = coord.idx * T::MERKLE_ARITY;
        let left_child_height = 0;
        let left_coord = Coord { height: left_child_height, idx: left_child_idx };

        let right_child_idx = left_child_idx + 1;
        let right_child_height = 0;
        let right_coord = Coord { height: right_child_height, idx: right_child_idx };

        if (!self.present_node.contains(&left_coord)) | (!self.present_node.contains(&right_coord)) {
            if self.contains_key_in_cache(coord) {
                BigMerkleTree::remove_node_from_cache(self, coord);
                if self.print_verbose { println!("remove_subtree_from_cache: Remove node from cache {:?}", coord); }
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
        if !self.present_node.contains(&coord) {
            if self.print_verbose { println!("node: Coord {:?} is an empty node.", coord); }
            return T::EMPTY_HASH_CST[coord.height];
        }
        if self.print_verbose { println!("node: Coord {:?} is not an empty node, retrieving node.", coord); }
        let res = self.get_from_cache(coord);

        // if not in the cache, recompute it
        if res == None {
            if self.print_verbose { println!("node: Node not in the cache, recomputing it"); }
            /* the node is not in the cache, compute it */
            let node_hash;
            if coord.height == 1 {
                /* get leaves to compute */
                let left_child_idx = coord.idx * T::MERKLE_ARITY;
                let left_child = self.get_from_db(left_child_idx);
                let left_hash: F;
                if let Some(i) = left_child {
                    left_hash = i;
                    if self.print_verbose { println!("node: Left leaf {} retrieved from db", left_child_idx); }
                } else {
                    left_hash = T::EMPTY_HASH_CST[0];
                    if self.print_verbose { println!("node: Left leaf {} is empty", left_child_idx); }
                }

                let right_child_idx = left_child_idx + 1;
                let right_child = self.get_from_db(right_child_idx);
                let right_hash: F;
                if let Some(i) = right_child {
                    right_hash = i;
                    if self.print_verbose { println!("node: Right leaf {} retrieved from db", right_child_idx); }
                } else {
                    right_hash = T::EMPTY_HASH_CST[0];
                    if self.print_verbose { println!("node: Right leaf {} is empty", right_child_idx); }
                }
                node_hash = merkle_tree::field_based_mht::smt::BigMerkleTree::<F, T, P>::poseidon_hash(left_hash, right_hash);
                if self.print_verbose { println!("node: hash of node {:?} computed. ", coord) };
            } else {
                let height_child = coord.height - 1;
                let left_child_idx = coord.idx * T::MERKLE_ARITY;
                let coord_left = Coord { height: height_child, idx: left_child_idx };
                let left_child_hash = BigMerkleTree::node(self, coord_left);

                let right_child_idx = left_child_idx + 1;
                let coord_right = Coord { height: height_child, idx: right_child_idx };
                let right_child_hash = BigMerkleTree::node(self, coord_right);

                node_hash = merkle_tree::field_based_mht::smt::BigMerkleTree::<F, T, P>::poseidon_hash(left_child_hash, right_child_hash);
                if self.print_verbose { println!("node: hash of node {:?} computed. ", coord) };
            }
            return node_hash;
        }

        res.unwrap()
    }

    pub fn insert_node_in_cache(&mut self, coord: Coord, hash: F) {
        self.insert_to_cache(coord, hash);
        //self.cache.insert(coord, hash);
    }

    pub fn remove_node_from_cache(&mut self, coord: Coord) {
        self.remove_from_cache(coord);
        //self.cache.remove(&coord);
    }

    pub fn poseidon_hash(x: F, y: F) -> F {
        let mut input = Vec::new();
        input.push(x);
        input.push(y);
        let hash = PoseidonHash::<F, P>::evaluate(&input[..]);
        hash.unwrap()
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
                self.present_node.remove(&coord);
            } else {
                self.insert_to_db(idx, hash.unwrap());
                self.present_node.insert(coord);
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
        let output_vec = merkle_tree::field_based_mht::smt::BigMerkleTree::<F, T, P>::batch_poseidon_hash(input_vec);
        // Place the computed hash in a cache_parallel
        let mut index_output_vec = 0;
        for coord in nodes_to_process_in_parallel[0].clone() {
            self.cache_parallel.insert(coord, output_vec[index_output_vec]);
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

                if self.cache_parallel.contains_key(&left_child_coord) {
                    left_hash = *self.cache_parallel.get(&left_child_coord).unwrap();
                } else {
                    left_hash = self.node(left_child_coord);
                }
                if self.cache_parallel.contains_key(&right_child_coord) {
                    right_hash = *self.cache_parallel.get(&right_child_coord).unwrap();
                } else {
                    right_hash = self.node(right_child_coord);
                }
                input_vec.push(left_hash);
                input_vec.push(right_hash);
                if self.present_node.contains(&left_child_coord) || self.present_node.contains(&right_child_coord){
                    self.present_node.insert(coord);
                } else {
                    self.present_node.remove(&coord);
                }
                if self.present_node.contains(&left_child_coord) && self.present_node.contains(&right_child_coord){
                    both_children_present.push(true);
                } else {
                    both_children_present.push(false);
                }
            }

            // Process the input_vec using batch Poseidon hash
            let output_vec = merkle_tree::field_based_mht::smt::BigMerkleTree::<F, T, P>::batch_poseidon_hash(input_vec);

            // Place the computed hash in a cache_parallel
            let mut index_output_vec = 0;
            for coord in nodes_to_process_in_parallel[height-1].clone() {
                self.cache_parallel.insert(coord, output_vec[index_output_vec]);
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

        self.root = *self.cache_parallel.get(&Coord{height:self.height,idx:0}).unwrap();
        self.root
    }

    pub fn process_leaves_normal(&mut self, lidx: Vec<OperationLeaf<F>>) -> F {

        for k in 0..lidx.len() {
            let x = lidx[k];
            let coord = x.coord;
            let action = x.action;
            let hash = x.hash;

            if action == Insert {
                self.insert_leaf(coord, hash.unwrap());
            } else {
                self.remove_leaf(coord);
            }
        }
        self.root
    }

}

pub type MNT4PoseidonSmt = BigMerkleTree<MNT4753Fr, MNT4753SmtPoseidonParameters, MNT4753PoseidonParameters>;

#[cfg(test)]
mod test {
    use crate::merkle_tree::field_based_mht::smt::{MNT4PoseidonSmt, MNT4PoseidonHash, OperationLeaf, Coord, ActionLeaf, SmtPoseidonParameters};
    use algebra::fields::mnt6753::Fr as MNT6753Fr;
    use algebra::fields::mnt4753::Fr as MNT4753Fr;
    use std::str::FromStr;
    use rand_xorshift::XorShiftRng;
    use rand::{SeedableRng, RngCore};
    use crate::merkle_tree::field_based_mht::{FieldBasedMerkleTreeConfig, FieldBasedMerkleHashTree};
    use std::time::Instant;
    use rand::rngs::OsRng;
    use crate::merkle_tree::field_based_mht::smt::parameters::MNT4753SmtPoseidonParameters;

    struct MNT4753FieldBasedMerkleTreeParams;

    impl FieldBasedMerkleTreeConfig for MNT4753FieldBasedMerkleTreeParams {
        const HEIGHT: usize = 6;
        type H = MNT4PoseidonHash;
    }

    type MNT4753FieldBasedMerkleTree = FieldBasedMerkleHashTree<MNT4753FieldBasedMerkleTreeParams>;

    #[test]
    fn process_leaves_mnt4_comp() {
        use algebra::{
            fields::mnt4753::Fr,
            UniformRand,
        };

        let num_leaves = 2usize.pow(23);
        //let num_leaves = 128;
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
            let mut smt1 = MNT4PoseidonSmt::new(num_leaves).unwrap();
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
            let mut smt2 = MNT4PoseidonSmt::new(num_leaves).unwrap();
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

        let mut smt = MNT4PoseidonSmt::new(num_leaves).unwrap();
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
        for _ in 30..32 {
            let f = Fr::zero();
            leaves.push(f);
        }
        let tree = MNT4753FieldBasedMerkleTree::new(&leaves).unwrap();

        println!("root_mt  = {:?}", tree.root.unwrap());
        println!("root_smt = {:?}", smt.root);
    }

    #[test]
    fn compare_merkle_trees_mnt4_1() {
        use algebra::{
            fields::mnt4753::Fr,
            UniformRand,
        };

        let mut rng = XorShiftRng::seed_from_u64(9174123u64);
        let num_leaves = 32;
        let mut leaves = Vec::new();
        for _ in 0..num_leaves {
            let f = Fr::rand(&mut rng);
            leaves.push(f);
        }
        let tree = MNT4753FieldBasedMerkleTree::new(&leaves).unwrap();

        let mut smt = MNT4PoseidonSmt::new(num_leaves).unwrap();
        let mut rng = XorShiftRng::seed_from_u64(9174123u64);
        for i in 0..num_leaves {
            let f = Fr::rand(&mut rng);
            smt.insert_leaf(
                Coord{height:0, idx:i},
                f,
            );
        }

        assert_eq!(tree.root(), smt.root, "Outputs of the Merkle trees for MNT4 do not match.");
    }

    #[test]
    fn compare_merkle_trees_mnt4_2() {
        use algebra::{
            fields::mnt4753::Fr, Field,
            //UniformRand,
        };

        let num_leaves = 32;
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
        for _ in 30..32 {
            let f = Fr::zero();
            leaves.push(f);
        }
        let tree = MNT4753FieldBasedMerkleTree::new(&leaves).unwrap();

        let mut smt = MNT4PoseidonSmt::new(num_leaves).unwrap();
        smt.insert_leaf(Coord{height:0, idx:0}, MNT4753Fr::from_str("1").unwrap());
        smt.insert_leaf(Coord{height:0, idx:9}, MNT4753Fr::from_str("2").unwrap());
        smt.insert_leaf(Coord{height:0, idx:16}, MNT4753Fr::from_str("10").unwrap());
        smt.insert_leaf(Coord{height:0, idx:29}, MNT4753Fr::from_str("3").unwrap());
        smt.remove_leaf(Coord{height:0, idx:16});

        assert_eq!(tree.root(), smt.root, "Outputs of the Merkle trees for MNT4 do not match.");
    }

}
