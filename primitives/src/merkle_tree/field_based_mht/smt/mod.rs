#![allow(dead_code)]

pub mod big_merkle_tree;
pub use self::big_merkle_tree::*;

pub mod big_lazy_merkle_tree;
pub use self::big_lazy_merkle_tree::*;

use algebra::{
    ToBytes, FromBytes, Field
};

use crate::merkle_tree::field_based_mht::FieldBasedMerkleTreeParameters;

use std::{
    io::{Write, Result as IoResult, Read},
    collections::{HashMap, HashSet},
    marker::PhantomData
};

#[derive(Debug, PartialEq, Eq, Hash, Copy, Clone)]
pub enum ActionLeaf {
    Insert,
    Remove,
}

#[derive(Debug, PartialEq, Eq, Hash, Copy, Clone)]
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

impl ToBytes for Coord {
    fn write<W: Write>(&self, mut writer: W) -> IoResult<()> {
        (self.height as u8).write(&mut writer)?;
        (self.idx as u32).write(&mut writer)
    }
}

impl FromBytes for Coord {
    fn read<R: Read>(mut reader: R) -> IoResult<Self> {
        let height = u8::read(&mut reader)? as usize;
        let idx = u32::read(&mut reader)? as usize;
        Ok(Self::new(height, idx))
    }
}

#[derive(Debug, PartialEq, Eq, Hash, Copy, Clone)]
// Action associated to the leaf
pub struct OperationLeaf <F: Field>{
    coord: Coord,
    action: ActionLeaf,
    hash: Option<F>,
}

impl<F: Field> OperationLeaf<F> {
    pub fn new(height: usize, idx: usize, action: ActionLeaf, hash: Option<F>) -> Self {
        Self { coord: Coord { height, idx }, action, hash}
    }
}

#[derive(Debug)]
pub(crate) struct BigMerkleTreeState<T: FieldBasedMerkleTreeParameters>{
    // the height of this tree
    height: usize,
    // stores the nodes of the path
    cache_path: HashMap<Coord, T::Data>,
    // indicates which nodes are present the Merkle tree
    present_node: HashSet<Coord>,
    // root of the Merkle tree
    root: T::Data,
    _parameters: PhantomData<T>
}

impl<T: FieldBasedMerkleTreeParameters> BigMerkleTreeState<T> {
    fn get_default_state(height: usize) -> Self {
        Self{
            height,
            cache_path: HashMap::new(),
            present_node: HashSet::new(),
            root: T::EMPTY_HASH_CST.unwrap().nodes[height],
            _parameters: PhantomData,
        }
    }
}

impl<T: FieldBasedMerkleTreeParameters> ToBytes for BigMerkleTreeState<T> {
    fn write<W: Write>(&self, mut writer: W) -> IoResult<()> {
        (self.height as u8).write(&mut writer)?;
        (self.cache_path.len() as u8).write(&mut writer)?;
        for (&coord, &fe) in self.cache_path.iter() {
            coord.write(&mut writer)?;
            fe.write(&mut writer)?;
        }

        (self.present_node.len() as u32).write(&mut writer)?;
        for &coord in self.present_node.iter() {
            coord.write(&mut writer)?;
        }

        self.root.write(&mut writer)
    }
}

impl<T: FieldBasedMerkleTreeParameters> FromBytes for BigMerkleTreeState<T> {
    fn read<R: Read>(mut reader: R) -> IoResult<Self> {
        let height = u8::read(&mut reader)? as usize;
        let cache_path_len = u8::read(&mut reader)? as usize;
        let mut cache_path = HashMap::new();
        for _ in 0..cache_path_len {
            let coord = Coord::read(&mut reader)?;
            let fe = T::Data::read(&mut reader)?;
            cache_path.insert(coord, fe);
        }

        let present_node_len = u32::read(&mut reader)? as usize;
        let mut present_node = HashSet::new();
        for _ in 0..present_node_len {
            let coord = Coord::read(&mut reader)?;
            present_node.insert(coord);
        }

        let root = T::Data::read(&mut reader)?;

        Ok(Self{height, cache_path, present_node, root, _parameters: PhantomData})
    }
}

