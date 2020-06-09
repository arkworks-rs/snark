use crate::merkle_tree::field_based_mht::smt::error::Error;
use std::collections::HashMap;
use std::collections::HashSet;

pub mod error;

pub const MERKLE_ARITY: usize = 2;
// Here it goes the precomputed hashes for the different levels of the Merkle tree
pub const EMPTY_HASH_CST: &'static[usize] = &[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15];

#[derive(Debug)]
pub struct BigMerkleTree {
    // the number of leaves
    width: usize,
    // the height of the Merkle tree
    height: usize,
    // stores the leaves
    db: HashMap<usize, usize>,
    // stores the cached nodes
    cache: HashMap<Coord, usize>,
    // stores the nodes of the path
    cache_path: HashMap<Coord, usize>,
    // indicates which nodes are not empty
    present_node: HashSet<Coord>,
    // root of the Merkle tree
    root: usize,
}

#[derive(Debug, PartialEq, Eq, Hash, Copy, Clone)]
// Coordinates system for identifying a node
pub struct Coord {
    height: usize,
    // height in the Merkle tree (0 -> leaves)
    idx: usize,             // the index of the node in that level
}

impl Coord {
    pub fn new(height: usize, idx: usize) -> Self {
        Self { height, idx }
    }
}

// Assumption: MERKLE_ARITY == 2
impl BigMerkleTree {
    pub fn new(width: usize) -> Result<Self, Error> {
        let height = width as f64;
        let height = height.log(MERKLE_ARITY as f64) as usize;
        let db = HashMap::new();
        let cache = HashMap::new();
        let cache_path = HashMap::new();
        let non_empty_list = HashSet::new();
        let root = EMPTY_HASH_CST[height];
        Ok(BigMerkleTree {
            width,
            height,
            db,
            cache,
            cache_path,
            present_node: non_empty_list,
            root,
        })
    }

    pub fn arity(&self) -> usize { MERKLE_ARITY }
    pub fn width(&self) -> usize { self.width }
    pub fn height(&self) -> usize { self.height }
    pub fn db(&self) -> &HashMap::<usize, usize> { &self.db }
    pub fn cache(&self) -> &HashMap::<Coord, usize> { &self.cache }

    // inserts a leaf in the Merkle tree
    // it updates the Merkle tree on the path from the leaf to the root
    pub fn insert_leaf(&mut self, idx: usize, leaf: usize) {
        let coord = Coord { height: 0, idx };
        // mark as non empty leaf
        self.present_node.insert(coord);
        println!("insert_leaf: insert leaf in non_empty set {:?}", coord);
        self.db.insert(idx, leaf);
        println!("insert_leaf: insert leaf in db {:?}", coord);
        self.cache_path.clear();
        self.cache_path.insert(coord, leaf);
        BigMerkleTree::update_tree(self, idx);
    }

    // removes a leaf in the Merkle tree
    // it updates the Merkle tree on the path from the leaf to the root
    pub fn remove_leaf(&mut self, idx: usize) {
        let coord = Coord { height: 0, idx };
        // take that leaf from the non-empty set
        self.present_node.remove(&coord);
        println!("remove_leaf: removes leaf from non_empty set {:?}", coord);
        self.cache_path.clear();
        // removes the leaf from the db
        let res = self.db.remove(&idx);
        println!("remove_leaf: removes leaf from db {:?}", coord);
        // if it was in the db, update the tree
        if res != None {
            BigMerkleTree::update_tree(self, idx);
        }
    }

    // it updates the tree visiting the parent nodes from the leaf to the root
    pub fn update_tree(&mut self, mut idx: usize) {

        // Process the node of level 1 with the inserted/removed leaf
        // check whether the hash corresponds to the left or right child
        let left_child_idx: usize;
        let right_child_idx: usize;
        let left_child: Option<&usize>;
        let right_child: Option<&usize>;
        let left_hash: usize;
        let right_hash: usize;
        let mut height = 0;

        if idx % MERKLE_ARITY == 0 {
            left_child_idx = idx;
            let coord = Coord { height, idx:left_child_idx };
            // get the left child from the cache_path
            let hash = self.cache_path.get(&coord);
            if let Some(i) = hash {
                left_hash = *i;
            } else {
                left_hash = EMPTY_HASH_CST[0];
            }
            right_child_idx = idx + 1;
            right_child = self.db.get(&right_child_idx);
            if let Some(i) = right_child {
                right_hash = *i;
            } else {
                right_hash = EMPTY_HASH_CST[0];
            }
        } else {
            right_child_idx = idx;
            let coord = Coord { height, idx: right_child_idx };
            // get the right child from the cache path
            let hash = self.cache_path.get(&coord);
            if let Some(i) = hash {
                right_hash = *i;
            } else {
                right_hash = EMPTY_HASH_CST[0];
            }
            left_child_idx = idx - 1;
            left_child = self.db.get(&left_child_idx);
            if let Some(i) = left_child {
                left_hash = *i;
            } else {
                left_hash = EMPTY_HASH_CST[0];
            }
        }

        height += 1;
        idx = idx / MERKLE_ARITY;
        let parent_coord = Coord { height, idx };

        let mut node_hash: usize;
        if (left_hash == EMPTY_HASH_CST[0]) & (right_hash == EMPTY_HASH_CST[0]) {
            node_hash = EMPTY_HASH_CST[height];
            println!("update_tree: processed leaves: left_hash (h: {}, idx: {}) = {}, right_hash (h: {}, idx: {}) = {}, node_hash = {}", 0, left_child_idx, left_hash, 0, right_child_idx, right_hash, node_hash);
            // insert the parent node into the cache_path
            self.cache_path.insert(parent_coord, node_hash);
            println!("-----> update_tree: visiting node {:?}", parent_coord);
            // Both children are empty leaves
            // remove the parent node from the presence set
            self.present_node.remove(&parent_coord.clone());
            println!("update_tree: remove node from empty_list {:?}", parent_coord);
            // remove node from cache
            BigMerkleTree::remove_node(self, parent_coord.clone());
            println!("update_tree: remove node from cache {:?}", parent_coord);
        } else {
            // compute the hash of the node with the hashes of the children
            node_hash = poseidon_hash(left_hash, right_hash);
            println!("update_tree: processed leaves: left_hash (h: {}, idx: {}) = {}, right_hash (h: {}, idx: {}) = {}, node_hash = {}", 0, left_child_idx, left_hash, 0, right_child_idx, right_hash, node_hash);
            // insert the parent node into the cache_path
            self.cache_path.insert(parent_coord, node_hash);
            println!("-----> update_tree: visiting node {:?}", parent_coord);
            // set the parent as present
            self.present_node.insert(parent_coord.clone());
            println!("update_tree: insert node in non_empty set {:?}", parent_coord);
            // Both children are not empty leaves
            if (left_hash != EMPTY_HASH_CST[0]) & (right_hash != EMPTY_HASH_CST[0]) {
                // cache the node
                BigMerkleTree::insert_node(self, parent_coord.clone(), node_hash);
                println!("update_tree: insert node in cache {:?}", parent_coord);
            }
        }

        // Process level >= 2
        while height != self.height {
            // go to parent node
            height += 1;
            let idx_child = idx;
            idx = idx / MERKLE_ARITY;
            let parent_coord = Coord { height, idx };
            println!("-----> update_tree: visiting node {:?}", parent_coord);

            // retrieve the left child
            let left_child_idx = parent_coord.idx * MERKLE_ARITY;
            let left_child_height = parent_coord.height - 1;
            let left_coord = Coord { height: left_child_height, idx: left_child_idx };
            println!("node coord left = {:?}", left_coord);
            let left_hash: usize;
            if left_child_idx == idx_child {
                // It was cached in the previous iteration. Get it from the cache_path
                left_hash = *self.cache_path.get(&left_coord).unwrap();
            } else {
                left_hash = BigMerkleTree::node(self, left_coord);
            }

            // retrieve the right child
            let right_child_idx = left_child_idx + 1;
            let right_child_height = left_child_height;
            let right_coord = Coord { height: right_child_height, idx: right_child_idx };
            println!("node coord right = {:?}", right_coord);
            let right_hash: usize;
            if right_child_idx == idx_child {
                // It was cached in the previous iteration. Get it from the cache_path
                right_hash = *self.cache_path.get(&right_coord).unwrap();
            } else {
                right_hash = BigMerkleTree::node(self, right_coord);
            }

            if (left_hash == EMPTY_HASH_CST[0]) & (right_hash == EMPTY_HASH_CST[0]) {
                node_hash = EMPTY_HASH_CST[height];
                // insert the parent node into the cache_path
                self.cache_path.insert(parent_coord, node_hash);
                println!("update_tree: processing nodes: left_hash (h: {}, idx: {}) = {}, right_hash (h: {}, idx: {}) = {}, node_hash = {}", left_child_height, left_child_idx, left_hash, right_child_height, right_child_idx, right_hash, node_hash);
                // remove node from non_empty set
                self.present_node.remove(&parent_coord.clone());
                println!("update_tree: remove node from empty_list {:?}", parent_coord);
            } else {
                // compute the hash of the parent node based on the hashes of the children
                node_hash = poseidon_hash(left_hash, right_hash);
                // insert the parent node into the cache_path
                self.cache_path.insert(parent_coord, node_hash);
                println!("update_tree: processing nodes: left_hash (h: {}, idx: {}) = {}, right_hash (h: {}, idx: {}) = {}, node_hash = {}", left_child_height, left_child_idx, left_hash, right_child_height, right_child_idx, right_hash, node_hash);
                if (left_hash != EMPTY_HASH_CST[left_child_height]) & (right_hash != EMPTY_HASH_CST[right_child_height]) {
                    // set the parent node as non_empty
                    self.present_node.insert(parent_coord.clone());
                    // children not empty leaves, then cache the parent node
                    BigMerkleTree::insert_node(self, parent_coord.clone(), node_hash);
                    println!("update_tree: insert node in cache {:?}", parent_coord);
                    // cache the children
                    BigMerkleTree::insert_node(self, left_coord.clone(), left_hash);
                    println!("update_tree: insert node in cache {:?}", left_coord);
                    BigMerkleTree::insert_node(self, right_coord.clone(), right_hash);
                    println!("update_tree: insert node in cache {:?}", right_coord);
                } else {
                    // either child not equal to empty leaves, include the parent node in a non-empty set
                    self.present_node.insert(parent_coord.clone());
                }
            }
            if (left_hash == EMPTY_HASH_CST[left_child_height]) | (right_hash == EMPTY_HASH_CST[right_child_height]) {
                if self.cache.contains_key(&parent_coord) {
                    // remove subtree from cache
                    BigMerkleTree::remove_subtree_from_cache(self, parent_coord);
                }
            }
        }
        self.root = node_hash;
    }

    pub fn remove_subtree_from_cache(&mut self, coord: Coord) {
        if coord.height == 1 {
            let left_child_idx = coord.idx * MERKLE_ARITY;
            let left_child_height = coord.height - 1;
            let left_coord = Coord { height: left_child_height, idx: left_child_idx };

            let right_child_idx = coord.idx * MERKLE_ARITY + 1;
            let right_child_height = coord.height - 1;
            let right_coord = Coord { height: right_child_height, idx: right_child_idx };

            if (!self.present_node.contains(&left_coord)) | (!self.present_node.contains(&right_coord)) {
                if self.cache.contains_key(&coord) {
                    BigMerkleTree::remove_node(self, coord);
                    println!("remove_subtree_from_cache: remove node from cache {:?}", coord);
                }
            }

            return;
        }
        if coord.height > 1 {
            let left_child_idx = coord.idx * MERKLE_ARITY;
            let left_child_height = coord.height - 1;
            let left_coord = Coord { height: left_child_height, idx: left_child_idx };

            let right_child_idx = coord.idx * MERKLE_ARITY + 1;
            let right_child_height = coord.height - 1;
            let right_coord = Coord { height: right_child_height, idx: right_child_idx };

            if (!self.present_node.contains(&left_coord)) | (!self.present_node.contains(&right_coord)) {
                if self.cache.contains_key(&coord) {
                    BigMerkleTree::remove_node(self, coord);
                    println!("remove_subtree_from_cache: remove node from cache {:?}", coord);
                }
            }

            let left_coord = Coord { height: coord.height - 1, idx: (coord.idx / MERKLE_ARITY) };
            let right_coord = Coord { height: coord.height - 1, idx: (coord.idx / MERKLE_ARITY + 1) };
            BigMerkleTree::remove_subtree_from_cache(self, left_coord);
            BigMerkleTree::remove_subtree_from_cache(self, right_coord);
        }
    }

    pub fn node(&mut self, coord: Coord) -> usize{
        //let coord = Coord{height, idx};
        // if the node is an empty node return the hash constant
        if !self.present_node.contains(&coord) {
            println!("Coord: {:?} is an empty node.", coord);
            return EMPTY_HASH_CST[coord.height];
        }
        println!("Coord: {:?} is NOT an empty node, retrieving node.", coord);
        let res = self.cache.get(&coord);

        // if not in the cache, recompute it
        if res == None {
            println!("Node not in the cache, recomputing it");
            /* the node is not in the cache, compute it */
            if coord.height == 1 {
                /* get leaves to compute */
                let left_child_idx = coord.idx * 2;
                let left_child = self.db.get(&left_child_idx);
                let left_hash:usize;
                if let Some(i) = left_child {
                    left_hash = *i;
                } else {
                    left_hash = 0;
                }

                let right_child_idx = left_child_idx + 1;
                let right_child = self.db.get(&right_child_idx);
                let right_hash:usize;
                if let Some(i) = right_child {
                    right_hash = *i;
                } else {
                    right_hash = 0;
                }
                let node_hash = poseidon_hash(left_hash, right_hash);
                println!("node: left_hash (h: {}, idx: {}) = {}, right_hash (h: {}, idx: {}) = {}, node_hash = {}", 0, left_child_idx, left_hash, 0, right_child_idx, right_hash, node_hash);

                return node_hash;
            } else {
                let height_child = coord.height - 1;
                let left_child_idx = coord.idx * 2;
                let coord_left = Coord{height:height_child, idx:left_child_idx};
                let left_child_hash = BigMerkleTree::node(self, coord_left);

                let right_child_idx = left_child_idx + 1;
                let coord_right = Coord{height:height_child, idx:right_child_idx};
                let right_child_hash = BigMerkleTree::node(self, coord_right);

                let node_hash = poseidon_hash(left_child_hash, right_child_hash);

                return node_hash;
            }
        }

        *res.unwrap()

    }

    pub fn node_is_empty(&self, height: usize, idx: usize) -> bool {
        /* if the node is not in the cache returns true */
        let coord = Coord{height, idx};
        !self.cache.contains_key(&coord)
    }

    pub fn insert_node(&mut self, coord: Coord, hash: usize) {
        self.cache.insert(coord, hash);
    }

    pub fn remove_node(&mut self, coord: Coord) {
        self.cache.remove(&coord);
    }

}

pub fn poseidon_hash (x: usize, y: usize) -> usize {
    if (x==0) & (y==0) {
        return 1;
    }
    x + y
}

#[cfg(test)]
mod test {
    use crate::merkle_tree::field_based_mht::smt::{BigMerkleTree};

    #[test]
    fn test_rocksdb() {
        let mut smt = BigMerkleTree::new(8).unwrap();
        println!("smt = {:?}", smt);

        smt.insert_leaf(
            0,
            10
        );

        println!("===============================================");

        smt.insert_leaf(
            2,
            20
        );

        println!("===============================================");

        smt.remove_leaf(
            2
        );


        println!("smt = {:?}", smt);
    }
}







// use rocksdb::{DB, Options};
// use std::sync::Arc;
// // NB: db is automatically closed at end of lifetime
//
// use algebra::fields::mnt4753::Fr as MNT4753Fr;
// use std::str::FromStr;
//
// use serde::{Serialize, Deserialize};
//
// #[derive(Debug)]
// pub struct BigMerkleTree {
//     width: usize,
//     height: usize,
//     db: Arc<DB>,
//     cache: Arc<DB>,
// }
//
// #[derive(Serialize, Deserialize, Debug)]
// struct field_elem {
//     x: MNT4753Fr,
// }
//
// pub fn test_rdb() {
//
//     let path = "/home/mkaihara/Documents/dev/ginger-lib/primitives/src/merkle_tree/field_based_mht/smt/rocksdb_storage";
//     {
//         let a = field_elem{x:MNT4753Fr::from_str("1").unwrap()};
//         let elem = bincode::serialize(&a).unwrap();
//
//         let db = DB::open_default(path).unwrap();
//         db.put(b"my key", elem).unwrap();
//         match db.get(b"my key") {
//             Ok(Some(value)) => {
//                 let e: field_elem = bincode::deserialize(&value).unwrap();
//                 println!("retrieved value {:?}", e);
//             },
//             Ok(None) => println!("value not found"),
//             Err(e) => println!("operational problem encountered: {}", e),
//         }
//         db.delete(b"my key").unwrap();
//     }
//     let _ = DB::destroy(&Options::default(), path);
// }
//
// pub fn insert(idx: usize, elem: MNT4753Fr) {
//
// }
//
// #[cfg(test)]
// mod test {
//     use crate::merkle_tree::field_based_mht::smt::test_rdb;
//
//     #[test]
//     fn test_rocksdb() {
//         test_rdb();
//     }
// }