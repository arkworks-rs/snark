#![allow(dead_code)]

pub mod big_merkle_tree;
pub mod big_lazy_merkle_tree;
pub mod error;

use crate::PoseidonHash;
use crate::crh::poseidon::parameters::{MNT4753PoseidonParameters, MNT6753PoseidonParameters};
use crate::merkle_tree::field_based_mht::FieldBasedMerkleTreeParameters;

use algebra::fields::mnt6753::Fr as MNT6753Fr;
use algebra::fields::mnt4753::Fr as MNT4753Fr;
use algebra::{PrimeField, ToBytes, FromBytes};

use serde::{Serialize,Deserialize};
use std::io::{Write, Result as IoResult, Read};
use std::collections::{HashMap, HashSet};
use std::marker::PhantomData;

pub type MNT4PoseidonHash = PoseidonHash<MNT4753Fr, MNT4753PoseidonParameters>;
pub type MNT6PoseidonHash = PoseidonHash<MNT6753Fr, MNT6753PoseidonParameters>;

#[derive(Debug, PartialEq, Eq, Hash, Copy, Clone)]
pub enum ActionLeaf {
    Insert,
    Remove,
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

impl ToBytes for Coord {
    fn write<W: Write>(&self, mut writer: W) -> IoResult<()> {
        (self.height as u8).write(&mut writer)?;
        (self.idx as u64).write(&mut writer)
    }
}

impl FromBytes for Coord {
    fn read<R: Read>(mut reader: R) -> IoResult<Self> {
        let height = u8::read(&mut reader)? as usize;
        let idx = u64::read(&mut reader)? as usize;
        Ok(Self::new(height, idx))
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

#[derive(Debug)]
pub(crate) struct BigMerkleTreeState<F: PrimeField, T: FieldBasedMerkleTreeParameters<Data = F>>{
    // the number of leaves
    width: usize,
    // stores the nodes of the path
    cache_path: HashMap<Coord, F>,
    // indicates which nodes are present the Merkle tree
    present_node: HashSet<Coord>,
    // root of the Merkle tree
    root: F,

    _parameters: PhantomData<T>
}

impl<F: PrimeField, T: FieldBasedMerkleTreeParameters<Data = F>> BigMerkleTreeState<F, T> {
    fn get_default_state(width: usize, height: usize) -> Self {
        Self{
            width,
            cache_path: HashMap::new(),
            present_node: HashSet::new(),
            root: T::EMPTY_HASH_CST[height],
            _parameters: PhantomData,
        }
    }
}

impl<F: PrimeField, T: FieldBasedMerkleTreeParameters<Data = F>> ToBytes for BigMerkleTreeState<F, T> {
    fn write<W: Write>(&self, mut writer: W) -> IoResult<()> {
        (self.width as u64).write(&mut writer)?;

        (self.cache_path.len() as u64).write(&mut writer)?;
        for (&coord, &fe) in self.cache_path.iter() {
            coord.write(&mut writer)?;
            fe.write(&mut writer)?;
        }

        (self.present_node.len() as u64).write(&mut writer)?;
        for &coord in self.present_node.iter() {
            coord.write(&mut writer)?;
        }

        self.root.write(&mut writer)
    }
}

impl<F: PrimeField, T: FieldBasedMerkleTreeParameters<Data = F>> FromBytes for BigMerkleTreeState<F, T> {
    fn read<R: Read>(mut reader: R) -> IoResult<Self> {
        let width = u64::read(&mut reader)? as usize;

        let cache_path_len = u64::read(&mut reader)? as usize;
        let mut cache_path = HashMap::new();
        for _ in 0..cache_path_len {
            let coord = Coord::read(&mut reader)?;
            let fe = F::read(&mut reader)?;
            cache_path.insert(coord, fe);
        }

        let present_node_len = u64::read(&mut reader)? as usize;
        let mut present_node = HashSet::new();
        for _ in 0..present_node_len {
            let coord = Coord::read(&mut reader)?;
            present_node.insert(coord);
        }

        let root = F::read(&mut reader)?;

        Ok(Self{width, cache_path, present_node, root, _parameters: PhantomData})
    }
}

