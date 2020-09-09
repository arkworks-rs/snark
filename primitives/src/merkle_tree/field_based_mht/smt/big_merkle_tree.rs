/* Sparse Big Merkle tree */
use algebra::{ToBytes, to_bytes, FromBytes};
use crate::{
    crh::{FieldBasedHash, FieldBasedHashParameters},
    merkle_tree::{
        field_based_mht::{
            FieldBasedMerkleTreeParameters,
            FieldBasedMerkleTreePath, FieldBasedBinaryMHTPath,
            smt::{
                Coord, OperationLeaf, ActionLeaf::Insert, BigMerkleTreeState,
            },
            check_precomputed_parameters
        }
    },
};

use rocksdb::{DB, Options};

use std::{
    marker::PhantomData, fs, path::Path, io::{Error, ErrorKind},
};

#[derive(Debug)]
pub struct BigMerkleTree<T: FieldBasedMerkleTreeParameters>{
    // if unset, all DBs and tree internal state will be deleted when an instance of this struct
    // gets dropped
    persistent: bool,
    // path to state required to restore the tree
    state_path: Option<String>,
    // tree in-memory state
    state: BigMerkleTreeState<T>,
    // the number of leaves
    width: usize,
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

impl<T: FieldBasedMerkleTreeParameters> Drop for BigMerkleTree<T> {
    fn drop(&mut self) {
        self.close();
    }
}

impl<T: FieldBasedMerkleTreeParameters> BigMerkleTree<T> {
    // Creates a new tree of specified `height`.
    // If `persistent` is specified, then DBs will be kept on disk and the tree state will be saved
    // so that the tree can be restored any moment later. Otherwise, no state will be saved on file
    // and the DBs will be deleted.
    pub fn new(
        height: usize,
        persistent: bool,
        state_path: Option<String>,
        path_db: String,
        path_cache: String
    ) -> Result<Self, Error> {
        let height = height - 1;
        assert!(check_precomputed_parameters::<T>(height));
        
        let rate = <<T::H  as FieldBasedHash>::Parameters as FieldBasedHashParameters>::R;

        assert_eq!(T::MERKLE_ARITY, 2); // For now we support only arity 2
        // Rate may also be smaller than the arity actually, but this assertion
        // is reasonable and simplify the design.
        assert_eq!(rate, T::MERKLE_ARITY);

        // If the tree must be persistent, that a path to which save the tree state must be
        // specified.
        if !persistent { assert!(state_path.is_none()) } else { assert!(state_path.is_some()) }
        
        let state = BigMerkleTreeState::<T>::get_default_state(height);
        let width = T::MERKLE_ARITY.pow(height as u32);

        let path_db = path_db;
        let database = DB::open_default(path_db.clone())
            .map_err(|e| Error::new(ErrorKind::Other, e))?;

        let path_cache = path_cache;
        let db_cache = DB::open_default(path_cache.clone())
            .map_err(|e| Error::new(ErrorKind::Other, e))?;

        Ok(BigMerkleTree {
            persistent,
            state_path,
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
    // in `new` function.
    pub fn load(
        persistent: bool,
        state_path: String,
        path_db: String,
        path_cache: String
    ) -> Result<Self, Error> {
        let rate = <<T::H  as FieldBasedHash>::Parameters as FieldBasedHashParameters>::R;
        assert_eq!(T::MERKLE_ARITY, 2); // For now we support only arity 2
        // Rate may also be smaller than the arity actually, but this assertion
        // is reasonable and simplify the design.
        assert_eq!(rate, T::MERKLE_ARITY);

        let state = {
            let state_file = fs::File::open(state_path.clone())?;
            BigMerkleTreeState::<T>::read(state_file)?
        };
        assert!(check_precomputed_parameters::<T>(state.height));
        let width = T::MERKLE_ARITY.pow(state.height as u32);

        let opening_options = Options::default();

        let path_db = path_db;
        let database = DB::open(&opening_options, path_db.clone())
            .map_err(|e| Error::new(ErrorKind::Other, e))?;

        let path_cache = path_cache;
        let db_cache = DB::open(&opening_options,path_cache.clone())
            .map_err(|e| Error::new(ErrorKind::Other, e))?;

        Ok(BigMerkleTree {
            persistent,
            state_path: Some(state_path),
            state,
            width,
            path_db,
            database,
            path_cache,
            db_cache,
            _parameters: PhantomData,
        })
    }

    pub fn close(&mut self) {
        if !self.persistent {

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
    
    pub fn set_persistency(&mut self, persistency: bool) {
        self.persistent = persistency;
    }

    /* ===============================================================================*/
    /* Cache operations */
    /*===============================================================================*/

    fn insert_to_cache(&self, coord: Coord, data: T::Data) {
        // Inserts the node into the cache

        let elem = to_bytes!(data).unwrap();
        let index = to_bytes!(coord).unwrap();
        self.db_cache.put(index, elem).unwrap();
    }

    fn get_from_cache(&self, coord:Coord) -> Option<T::Data> {
        // Retrieves the node from the cache
        // If the node is in the cache, it returns the hash as an option
        // If the node is not in the cache, or there is an error, it returns none

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
        // Deletes the node from the cache
        // If the node was in the cache, it deletes the node and returns the hash as an option
        // If the node was not present in the cache, or if there was an error, it returns none

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

    /* ===============================================================================*/
    /* Leaves operations */
    /*===============================================================================*/

    fn insert_to_db(&self, idx: usize, data: T::Data) {
        // Inserts the leaf to the db

        let elem = to_bytes!(data).unwrap();
        let index = to_bytes!(idx as u32).unwrap();
        self.database.put(index, elem).unwrap();
    }

    fn get_from_db(&self, idx: usize) -> Option<T::Data>{
        // Retrieves the leaf from the db

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
        // Deletes the leaf from the db
        // If the leaf was in the db, it return the hash as an option
        // If the leaf was not present, or if there is an error, it returns none

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

    /* ===============================================================================*/
    /* Merkle tree operations */
    /*===============================================================================*/

    pub fn insert_leaf(&mut self, coord: Coord, leaf: T::Data) {
        // Inserts a leaf in the Merkle tree.
        // Updates the Merkle tree on the path from the leaf to the root

        // check that the index of the leaf to be inserted is less than the width of the Merkle tree
        println!("idx: {}", coord.idx);
        println!("width: {}", self.width);
        assert!(coord.idx < self.width, "Leaf index out of bound.");
        // check that the coordinates of the node corresponds to the leaf level
        assert_eq!(coord.height, 0, "Coord of the node does not correspond to leaf level");

        if self.state.present_node.contains(&coord) {
            let old_hash = self.get_from_db(coord.idx).unwrap();
            if old_hash != leaf {
                self.insert_to_db(coord.idx, leaf);
                self.state.cache_path.clear();
                self.state.cache_path.insert(coord, leaf);
                BigMerkleTree::update_tree(self, coord);
            }
        } else {
            // mark as non empty leaf
            self.state.present_node.insert(coord);
            self.insert_to_db(coord.idx, leaf);
            self.state.cache_path.clear();
            self.state.cache_path.insert(coord, leaf);
            BigMerkleTree::update_tree(self, coord);
        }
    }

    // Removes a leaf in the Merkle tree.
    // Updates the Merkle tree on the path from the leaf to the root
    pub fn remove_leaf(&mut self, coord: Coord) {

        // check that the index of the leaf to be inserted is less than the width of the Merkle tree
        assert!(coord.idx < self.width, "Leaf index out of bound.");
        // check that the coordinates of the node corresponds to the leaf level
        assert_eq!(coord.height, 0, "Coord of the node does not correspond to leaf level");

        // take that leaf from the non-empty set
        self.state.present_node.remove(&coord);
        self.state.cache_path.clear();
        self.state.cache_path.insert(coord, T::EMPTY_HASH_CST.unwrap().nodes[0]);
        // removes the leaf from the db
        let res = self.remove_from_db(coord.idx);
        // if it was in the db, update the tree
        if res != None {
            BigMerkleTree::update_tree(self, coord);
        }
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
            let sibling = if self.state.present_node.contains(&sibling_coord) {
                // If it's not empty
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

            // go up one level
            height += 1;
            node_idx = node_idx / T::MERKLE_ARITY; // compute the index of the parent
        }
        assert_eq!(self.node(Coord { height, idx: node_idx }), self.state.root);

        return FieldBasedBinaryMHTPath::<T>::new(path);
    }

    pub fn is_leaf_empty(&self, coord: Coord) -> bool {
        // check that the index of the leaf is less than the width of the Merkle tree
        assert!(coord.idx < self.width, "Leaf index out of bound.");
        // check that the coordinates of the node corresponds to the leaf level
        assert_eq!(coord.height, 0, "Coord of the node does not correspond to leaf level");

        !self.state.present_node.contains(&coord)
    }

    // Updates the tree visiting the parent nodes from the leaf to the root
    // Calculates the hash and caches it
    fn update_tree(&mut self, coord: Coord) {

        // Process the node of level 1 with the inserted/removed leaf
        // check whether the hash corresponds to the left or right child
        let mut idx = coord.idx;
        let mut height = 0;
        let empty_nodes = T::EMPTY_HASH_CST.unwrap().nodes;
        let empty_node = empty_nodes[height];

        assert_eq!(T::MERKLE_ARITY,2, "Arity of the Merkle tree is not 2.");

        // Get left and right child at leaves level
        let (left_child_idx, right_child_idx, is_coord_left) = if idx % T::MERKLE_ARITY == 0 {
            (idx, idx + 1, true)
        } else {
            (idx - 1, idx, false)
        };

        let left_child_coord = Coord { height, idx: left_child_idx };
        let right_child_coord = Coord { height, idx: right_child_idx };

        let (left_hash, right_hash) = if is_coord_left {
            (
                *self.state.cache_path.get(&left_child_coord).unwrap_or(&empty_node),
                self.get_from_db(right_child_idx).unwrap_or(empty_node)
            )
        } else {
            (
                self.get_from_db(left_child_idx).unwrap_or(empty_node),
                *self.state.cache_path.get(&right_child_coord).unwrap_or(&empty_node)
            )
        };

        // go up one level
        height += 1;
        idx = idx / T::MERKLE_ARITY;
        let parent_coord = Coord { height, idx };

        let mut node_hash: T::Data;
        if (!self.state.present_node.contains(&left_child_coord)) & (!self.state.present_node.contains(&right_child_coord)) {
            // if both children are empty

            node_hash = empty_nodes[height];

            // insert the parent node into the cache_path
            self.state.cache_path.insert(parent_coord, node_hash);

            // Both children are empty leaves
            // remove the parent node from the presence set
            self.state.present_node.remove(&parent_coord.clone());
            // remove node from cache
            self.remove_from_cache(parent_coord.clone());
        } else {

            // compute the hash of the node with the hashes of the children
            node_hash = Self::field_hash(&left_hash, &right_hash);
            // insert the parent node into the cache_path
            self.state.cache_path.insert(parent_coord, node_hash);
            // set the parent as present
            self.state.present_node.insert(parent_coord.clone());
            // Both children are not empty leaves
            if (self.state.present_node.contains(&left_child_coord)) & (self.state.present_node.contains(&right_child_coord)) {
                // cache the node
                self.insert_to_cache(parent_coord.clone(), node_hash);
            }
        }

        // Process level >= 2
        while height != self.state.height {

            let left_child_height = height;
            let right_child_height = height;

            // go to parent node
            height += 1;
            let idx_child = idx;
            idx = idx / T::MERKLE_ARITY;
            let parent_coord = Coord { height, idx };

            // retrieve the children
            let left_child_idx = parent_coord.idx * T::MERKLE_ARITY;
            let right_child_idx = left_child_idx + 1;

            let left_child_coord = Coord { height: left_child_height, idx: left_child_idx };
            let right_child_coord = Coord { height: right_child_height, idx: right_child_idx };

            let (left_hash, right_hash) = if left_child_idx == idx_child {
                (
                    *self.state.cache_path.get(&left_child_coord).unwrap(),
                    BigMerkleTree::node(self, right_child_coord)
                )
            } else {
                (
                    BigMerkleTree::node(self, left_child_coord),
                    *self.state.cache_path.get(&right_child_coord).unwrap()
                )
            };

            if (!self.state.present_node.contains(&left_child_coord)) & (!self.state.present_node.contains(&right_child_coord)) {
                // both children are empty => parent as well

                node_hash = empty_nodes[height];
                // insert the parent node into the cache_path
                self.state.cache_path.insert(parent_coord, node_hash);
                // remove node from non_empty set
                self.state.present_node.remove(&parent_coord.clone());
            } else {
                // at least one is non-empty

                // compute the hash of the parent node based on the hashes of the children
                node_hash = Self::field_hash(&left_hash, &right_hash);
                // insert the parent node into the cache_path
                self.state.cache_path.insert(parent_coord, node_hash);

                if self.state.present_node.contains(&left_child_coord) & self.state.present_node.contains(&right_child_coord) {
                    // both children are present

                    // set the parent node as non_empty
                    self.state.present_node.insert(parent_coord.clone());
                    // children not empty leaves, then cache the parent node
                    self.insert_to_cache(parent_coord.clone(), node_hash);
                    // cache the children
                    self.insert_to_cache(left_child_coord.clone(), left_hash);
                    self.insert_to_cache(right_child_coord.clone(), right_hash);
                } else {
                    // one is empty and the other is non-empty child, include the parent node in a non-empty set
                    self.state.present_node.insert(parent_coord.clone());
                }
            }
            if (!self.state.present_node.contains(&left_child_coord)) | (!self.state.present_node.contains(&right_child_coord)) {
                // at least one child is empty

                // remove subtree from cache
                BigMerkleTree::remove_subtree_from_cache(self, parent_coord);

            }
        }
        self.state.root = node_hash;
    }

    pub fn get_root(&self) -> T::Data {
        self.state.root.clone()
    }

    pub fn height(&self) -> usize { self.state.height + 1 }

    fn remove_subtree_from_cache(&mut self, coord: Coord) {

        assert_eq!(T::MERKLE_ARITY,2, "Arity of the Merkle tree is not 2.");

        assert!(coord.height > 1);

        // remove the parent node from the cache
        self.remove_from_cache(coord.clone());

        let left_child_idx = coord.idx * T::MERKLE_ARITY;
        let left_child_height = coord.height - 1;
        let left_coord = Coord { height: left_child_height, idx: left_child_idx };

        if self.state.present_node.contains(&left_coord) {
            // go to the left child
            let ll_child_idx = left_child_idx * T::MERKLE_ARITY;
            let ll_child_height = left_child_height - 1;
            let ll_coord = Coord { height: ll_child_height, idx: ll_child_idx };

            let lr_child_idx = ll_child_idx + 1;
            let lr_child_height = ll_child_height;
            let lr_coord = Coord { height: lr_child_height, idx: lr_child_idx };

            if (!self.state.present_node.contains(&ll_coord)) | (!self.state.present_node.contains(&lr_coord)) {
                self.remove_from_cache(left_coord);
            }
        }

        let right_child_idx = left_child_idx + 1;
        let right_child_height = left_child_height;
        let right_coord = Coord { height: right_child_height, idx: right_child_idx };

        if self.state.present_node.contains(&right_coord) {
            // go to the right child

            let rl_child_idx = right_child_idx * T::MERKLE_ARITY;
            let rl_child_height = right_child_height - 1;
            let rl_coord = Coord { height: rl_child_height, idx: rl_child_idx };

            let rr_child_idx = rl_child_idx + 1;
            let rr_child_height = rl_child_height;
            let rr_coord = Coord { height: rr_child_height, idx: rr_child_idx };

            if (!self.state.present_node.contains(&rl_coord)) | (!self.state.present_node.contains(&rr_coord)) {
                self.remove_from_cache(right_coord);
            }
        }
        return;
    }

    fn node(&mut self, coord: Coord) -> T::Data {
        // Returns the hash value associated to the node.
        // If the node is in the cache, it retrieves from it.
        // If not, recomputes it.
        // Only used for nodes of level >= 1 (not for leaves).


        assert_eq!(T::MERKLE_ARITY,2, "Arity of the Merkle tree is not 2.");

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
                let left_child_hash = BigMerkleTree::node(self, coord_left);

                let right_child_idx = left_child_idx + 1;
                let coord_right = Coord { height: height_child, idx: right_child_idx };
                let right_child_hash = BigMerkleTree::node(self, coord_right);

                node_hash = Self::field_hash(&left_child_hash, &right_child_hash);
            }
            return node_hash;
        }

        res.unwrap()
    }

    fn field_hash(x: &T::Data, y: &T::Data) -> T::Data {
        <T::H as FieldBasedHash>::init(None)
            .update(x.clone())
            .update(y.clone())
            .finalize()
    }

    pub fn process_leaves_normal(&mut self, lidx: Vec<OperationLeaf<T::Data>>) -> T::Data {

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
        self.state.root.clone()
    }
}


#[cfg(test)]
mod test {

    use algebra::{
        biginteger::BigInteger768,
        fields::{mnt6753::Fr as MNT6753Fr, mnt4753::Fr as MNT4753Fr},
        Field, UniformRand
    };

    use crate::{
        crh::{
            MNT4PoseidonHash, MNT6PoseidonHash
        },
        merkle_tree::field_based_mht::{
            naive:: NaiveMerkleTree,
            smt::{OperationLeaf, Coord, ActionLeaf, BigMerkleTree},
            poseidon::{
                MNT4753_MHT_POSEIDON_PARAMETERS, MNT6753_MHT_POSEIDON_PARAMETERS
            },
            FieldBasedMerkleTreeParameters, FieldBasedMerkleTreePrecomputedEmptyConstants,
            FieldBasedMerkleTreePath
        },

    };

    use rand_xorshift::XorShiftRng;
    use rand::SeedableRng;

    use std::{
        str::FromStr, path::Path
    };

    #[derive(Clone, Debug)]
    struct MNT4753FieldBasedMerkleTreeParams;
    impl FieldBasedMerkleTreeParameters for MNT4753FieldBasedMerkleTreeParams {
        type Data = MNT4753Fr;
        type H = MNT4PoseidonHash;
        const MERKLE_ARITY: usize = 2;
        const EMPTY_HASH_CST: Option<FieldBasedMerkleTreePrecomputedEmptyConstants<'static, Self::H>> =
            Some(MNT4753_MHT_POSEIDON_PARAMETERS);
    }

    type MNT4753FieldBasedMerkleTree = NaiveMerkleTree<MNT4753FieldBasedMerkleTreeParams>;
    type MNT4PoseidonSMT = BigMerkleTree<MNT4753FieldBasedMerkleTreeParams>;

    #[derive(Clone, Debug)]
    struct MNT6753FieldBasedMerkleTreeParams;
    impl FieldBasedMerkleTreeParameters for MNT6753FieldBasedMerkleTreeParams {
        type Data = MNT6753Fr;
        type H = MNT6PoseidonHash;
        const MERKLE_ARITY: usize = 2;
        const EMPTY_HASH_CST: Option<FieldBasedMerkleTreePrecomputedEmptyConstants<'static, Self::H>> =
            Some(MNT6753_MHT_POSEIDON_PARAMETERS);
    }

    type MNT6753FieldBasedMerkleTree = NaiveMerkleTree<MNT6753FieldBasedMerkleTreeParams>;
    type MNT6PoseidonSMT = BigMerkleTree<MNT6753FieldBasedMerkleTreeParams>;

    const TEST_HEIGHT: usize = 6;

    #[test]
    fn compare_merkle_trees_mnt4_1() {

        let mut leaves_to_process: Vec<OperationLeaf<MNT4753Fr>> = Vec::new();

        leaves_to_process.push(OperationLeaf { coord: Coord { height: 0, idx: 0 }, action: ActionLeaf::Insert, hash: Some(MNT4753Fr::from_str("1").unwrap()) });
        leaves_to_process.push(OperationLeaf { coord: Coord { height: 0, idx: 9 }, action: ActionLeaf::Insert, hash: Some(MNT4753Fr::from_str("2").unwrap()) });
        leaves_to_process.push(OperationLeaf { coord: Coord { height: 0, idx: 16 }, action: ActionLeaf::Remove, hash: Some(MNT4753Fr::from_str("3").unwrap()) });
        leaves_to_process.push(OperationLeaf { coord: Coord { height: 0, idx: 29 }, action: ActionLeaf::Insert, hash: Some(MNT4753Fr::from_str("3").unwrap()) });
        leaves_to_process.push(OperationLeaf { coord: Coord { height: 0, idx: 16 }, action: ActionLeaf::Remove, hash: Some(MNT4753Fr::from_str("3").unwrap()) });

        let mut smt = MNT4PoseidonSMT::new(
            TEST_HEIGHT,
            false,
            None,
            String::from("./db_leaves_compare_merkle_trees_mnt4_1"),
            String::from("./db_cache_compare_merkle_trees_mnt4_1")
        ).unwrap();
        smt.process_leaves_normal(leaves_to_process);

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
        for _ in 30..32 {
            let f = MNT4753Fr::zero();
            leaves.push(f);
        }
        let mut tree = MNT4753FieldBasedMerkleTree::new(TEST_HEIGHT);
        tree.append(&leaves).unwrap();

        assert_eq!(tree.root(), smt.state.root, "Outputs of the Merkle trees for MNT4 do not match.");
    }

    #[test]
    fn compare_merkle_trees_mnt4_2() {

        let mut rng = XorShiftRng::seed_from_u64(9174123u64);
        let num_leaves = 32;
        let mut leaves = Vec::new();
        for _ in 0..num_leaves {
            let f = MNT4753Fr::rand(&mut rng);
            leaves.push(f);
        }
        let mut tree = MNT4753FieldBasedMerkleTree::new(TEST_HEIGHT);
        tree.append(&leaves).unwrap();

        let mut smt = MNT4PoseidonSMT::new(
            TEST_HEIGHT,
            false,
            None,
            String::from("./db_leaves_compare_merkle_trees_mnt4_2"),
            String::from("./db_cache_compare_merkle_trees_mnt4_2")
        ).unwrap();
        let mut rng = XorShiftRng::seed_from_u64(9174123u64);
        for i in 0..num_leaves {
            let f = MNT4753Fr::rand(&mut rng);
            smt.insert_leaf(
                Coord{height:0, idx:i},
                f,
            );
        }

        assert_eq!(tree.root(), smt.state.root, "Outputs of the Merkle trees for MNT4 do not match.");
    }

    #[test]
    fn compare_merkle_trees_mnt4_3() {

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
        for _ in 30..32 {
            let f = MNT4753Fr::zero();
            leaves.push(f);
        }

        let mut tree = MNT4753FieldBasedMerkleTree::new(TEST_HEIGHT);
        tree.append(&leaves).unwrap();

        let mut smt = MNT4PoseidonSMT::new(
            TEST_HEIGHT,
            false,
            None,
            String::from("./db_leaves_compare_merkle_trees_mnt4_3"),
            String::from("./db_cache_compare_merkle_trees_mnt4_3")
        ).unwrap();
        
        smt.insert_leaf(Coord{height:0, idx:0}, MNT4753Fr::from_str("1").unwrap());
        smt.insert_leaf(Coord{height:0, idx:9}, MNT4753Fr::from_str("2").unwrap());
        smt.insert_leaf(Coord{height:0, idx:16}, MNT4753Fr::from_str("10").unwrap());
        smt.insert_leaf(Coord{height:0, idx:29}, MNT4753Fr::from_str("3").unwrap());
        smt.remove_leaf(Coord{height:0, idx:16});

        assert_eq!(tree.root(), smt.state.root, "Outputs of the Merkle trees for MNT4 do not match.");
    }

    #[test]
    fn compare_merkle_trees_mnt6_1() {

        let mut leaves_to_process: Vec<OperationLeaf<MNT6753Fr>> = Vec::new();

        leaves_to_process.push(OperationLeaf { coord: Coord { height: 0, idx: 0 }, action: ActionLeaf::Insert, hash: Some(MNT6753Fr::from_str("1").unwrap()) });
        leaves_to_process.push(OperationLeaf { coord: Coord { height: 0, idx: 9 }, action: ActionLeaf::Insert, hash: Some(MNT6753Fr::from_str("2").unwrap()) });
        leaves_to_process.push(OperationLeaf { coord: Coord { height: 0, idx: 16 }, action: ActionLeaf::Remove, hash: Some(MNT6753Fr::from_str("3").unwrap()) });
        leaves_to_process.push(OperationLeaf { coord: Coord { height: 0, idx: 29 }, action: ActionLeaf::Insert, hash: Some(MNT6753Fr::from_str("3").unwrap()) });
        leaves_to_process.push(OperationLeaf { coord: Coord { height: 0, idx: 16 }, action: ActionLeaf::Remove, hash: Some(MNT6753Fr::from_str("3").unwrap()) });

        let mut smt = MNT6PoseidonSMT::new(
            TEST_HEIGHT,
            false,
            None,
            String::from("./db_leaves_compare_merkle_trees_mnt6_1"),
            String::from("./db_cache_compare_merkle_trees_mnt6_1")
        ).unwrap();
        smt.process_leaves_normal(leaves_to_process);

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
        for _ in 30..32 {
            let f = MNT6753Fr::zero();
            leaves.push(f);
        }
        let mut tree = MNT6753FieldBasedMerkleTree::new(TEST_HEIGHT);
        tree.append(&leaves).unwrap();

        assert_eq!(tree.root(), smt.state.root, "Outputs of the Merkle trees for MNT6 do not match.");
    }

    #[test]
    fn compare_merkle_trees_mnt6_2() {

        let mut rng = XorShiftRng::seed_from_u64(9174123u64);
        let num_leaves = 32;
        let mut leaves = Vec::new();
        for _ in 0..num_leaves {
            let f = MNT6753Fr::rand(&mut rng);
            leaves.push(f);
        }

        let mut tree = MNT6753FieldBasedMerkleTree::new(TEST_HEIGHT);
        tree.append(&leaves).unwrap();

        let mut smt = MNT6PoseidonSMT::new(
            TEST_HEIGHT,
            false,
            None,
            String::from("./db_leaves_compare_merkle_trees_mnt6_2"),
            String::from("./db_cache_compare_merkle_trees_mnt6_2")
        ).unwrap();
        let mut rng = XorShiftRng::seed_from_u64(9174123u64);
        for i in 0..num_leaves {
            let f = MNT6753Fr::rand(&mut rng);
            smt.insert_leaf(
                Coord{height:0, idx:i},
                f,
            );
        }

        assert_eq!(tree.root(), smt.state.root, "Outputs of the Merkle trees for MNT6 do not match.");
    }

    #[test]
    fn compare_merkle_trees_mnt6_3() {

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
        for _ in 30..32 {
            let f = MNT6753Fr::zero();
            leaves.push(f);
        }
        let mut tree = MNT6753FieldBasedMerkleTree::new(TEST_HEIGHT);
        tree.append(&leaves).unwrap();

        let mut smt = MNT6PoseidonSMT::new(
            TEST_HEIGHT,
            false,
            None,
            String::from("./db_leaves_compare_merkle_trees_mnt6_3"),
            String::from("./db_cache_compare_merkle_trees_mnt6_3")).unwrap();

        smt.insert_leaf(Coord{height:0, idx:0}, MNT6753Fr::from_str("1").unwrap());
        smt.insert_leaf(Coord{height:0, idx:9}, MNT6753Fr::from_str("2").unwrap());
        smt.insert_leaf(Coord{height:0, idx:16}, MNT6753Fr::from_str("10").unwrap());
        smt.insert_leaf(Coord{height:0, idx:29}, MNT6753Fr::from_str("3").unwrap());
        smt.remove_leaf(Coord{height:0, idx:16});

        println!("{:?}", smt.state.root);

        assert_eq!(tree.root(), smt.state.root, "Outputs of the Merkle trees for MNT6 do not match.");
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
            let mut smt = MNT4PoseidonSMT::new(
                TEST_HEIGHT,
                true,
                Some(String::from("./persistency_test_info")),
                String::from("./db_leaves_persistency_test_info"),
                String::from("./db_cache_persistency_test_info")
            ).unwrap();

            //Insert some leaves in the tree
            smt.insert_leaf(Coord{height:0, idx:0}, MNT4753Fr::from_str("1").unwrap());
            smt.insert_leaf(Coord{height:0, idx:9}, MNT4753Fr::from_str("2").unwrap());

            // smt gets dropped but its info should be saved
        }
        // files and directories should have been created
        assert!(Path::new("./persistency_test_info").exists());
        assert!(Path::new("./db_leaves_persistency_test_info").exists());
        assert!(Path::new("./db_cache_persistency_test_info").exists());

        // create a non-persistent smt in another scope by restoring the previous one
        {
            let mut smt = MNT4PoseidonSMT::load(
                false,
                String::from("./persistency_test_info"),
                String::from("./db_leaves_persistency_test_info"),
                String::from("./db_cache_persistency_test_info")
            ).unwrap();

            // insert other leaves
            smt.insert_leaf(Coord { height: 0, idx: 16 }, MNT4753Fr::from_str("10").unwrap());
            smt.insert_leaf(Coord { height: 0, idx: 29 }, MNT4753Fr::from_str("3").unwrap());

            // if truly state has been kept, then the equality below must pass, since `root` was
            // computed in one go with another smt
            assert_eq!(root, smt.state.root);

            // smt gets dropped and state and dbs are deleted
        }

        // files and directories should have been deleted
        assert!(!Path::new("./persistency_test_info").exists());
        assert!(!Path::new("./db_leaves_persistency_test_info").exists());
        assert!(!Path::new("./db_cache_persistency_test_info").exists());
    }

    #[test]
    fn merkle_tree_path_test_mnt4() {

        let num_leaves = 32;
        let mut leaves = Vec::with_capacity(num_leaves);
        let mut rng = XorShiftRng::seed_from_u64(1231275789u64);

        let mut smt = MNT4PoseidonSMT::new(
            TEST_HEIGHT,
            false,
            None,
            String::from("./db_leaves_merkle_tree_path_test_mnt4"),
            String::from("./db_cache_merkle_tree_path_test_mnt4")
        ).unwrap();

        // Generate random leaves, half of which empty
        for i in 0..num_leaves/2 {
            let leaf = MNT4753Fr::rand(&mut rng);
            leaves.push(leaf);
            smt.insert_leaf(Coord{height: 0, idx: i}, leaf);
        }

        for i in num_leaves/2..num_leaves {
            let leaf = MNT4753Fr::zero();
            leaves.push(leaf);
            smt.insert_leaf(Coord{height: 0, idx: i}, leaf);
        }

        // Compute the root of the tree, and do the same for a NaiveMHT, used here as reference
        let mut naive_tree = MNT4753FieldBasedMerkleTree::new(TEST_HEIGHT);
        naive_tree.append(&leaves).unwrap();
        let root = smt.get_root();
        let naive_root = naive_tree.root();
        assert_eq!(root, naive_root);

        for i in 0..num_leaves {

            // Create and verify a FieldBasedMHTPath
            let path = smt.get_merkle_path(Coord{height: 0, idx: i});
            assert!(path.verify( smt.height(), &leaves[i], &root).unwrap());

            // Create and verify a Naive path
            let naive_path = naive_tree.generate_proof(i, &leaves[i]).unwrap();
            assert!(naive_path.verify(naive_tree.height(), &leaves[i], &naive_root).unwrap());

            // Assert the two paths are equal
            assert_eq!(path, naive_path);
        }
    }

    #[test]
    fn merkle_tree_path_test_mnt6() {

        let num_leaves = 32;
        let mut leaves = Vec::with_capacity(num_leaves);
        let mut rng = XorShiftRng::seed_from_u64(1231275789u64);

        let mut smt = MNT6PoseidonSMT::new(
            TEST_HEIGHT,
            false,
            None,
            String::from("./db_leaves_merkle_tree_path_test_mnt6"),
            String::from("./db_cache_merkle_tree_path_test_mnt6")
        ).unwrap();

        // Generate random leaves, half of which empty
        for i in 0..num_leaves/2 {
            let leaf = MNT6753Fr::rand(&mut rng);
            leaves.push(leaf);
            smt.insert_leaf(Coord{height: 0, idx: i}, leaf);
        }

        for i in num_leaves/2..num_leaves {
            let leaf = MNT6753Fr::zero();
            leaves.push(leaf);
            smt.insert_leaf(Coord{height: 0, idx: i}, leaf);
        }

        // Compute the root of the tree, and do the same for a NaiveMHT, used here as reference
        let mut naive_tree = MNT6753FieldBasedMerkleTree::new(TEST_HEIGHT);
        naive_tree.append(&leaves).unwrap();
        let root = smt.get_root();
        let naive_root = naive_tree.root();
        assert_eq!(root, naive_root);

        for i in 0..num_leaves {

            // Create and verify a FieldBasedMHTPath
            let path = smt.get_merkle_path(Coord{height: 0, idx: i});
            assert!(path.verify( smt.height(), &leaves[i], &root).unwrap());

            // Create and verify a Naive path
            let naive_path = naive_tree.generate_proof(i, &leaves[i]).unwrap();
            assert!(naive_path.verify(naive_tree.height(), &leaves[i], &naive_root ).unwrap());

            // Assert the two paths are equal
            assert_eq!(path, naive_path);
        }
    }
}