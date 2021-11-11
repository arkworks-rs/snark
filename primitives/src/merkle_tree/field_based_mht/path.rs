use crate::{crh::*, field_based_mht::*};
use algebra::{serialize::*, SemanticallyValid};
use std::{
    clone::Clone,
    convert::TryFrom,
    io::{Read, Result as IoResult, Write},
};

/// An implementation of the FieldBasedMerkleTreePath trait, for a given FieldBasedHash and
/// FieldBasedMerkleTree with arbitrary arity.
/// TODO: Test for arity > 2
#[derive(Derivative)]
#[derivative(
    Clone(bound = ""),
    Debug(bound = ""),
    Default(bound = ""),
    Eq(bound = "")
)]
#[derive(Serialize, Deserialize, CanonicalSerialize, CanonicalDeserialize)]
pub struct FieldBasedMHTPath<T: FieldBasedMerkleTreeParameters> {
    path: Vec<(Vec<<T::H as FieldBasedHash>::Data>, usize)>,
}

impl<T: FieldBasedMerkleTreeParameters> SemanticallyValid for FieldBasedMHTPath<T> {
    fn is_valid(&self) -> bool {
        for (fes, pos) in self.path.iter() {
            if fes.len() != T::MERKLE_ARITY - 1 || pos >= &T::MERKLE_ARITY || !fes.is_valid() {
                return false;
            }
        }
        true
    }
}

impl<T: FieldBasedMerkleTreeParameters> PartialEq for FieldBasedMHTPath<T> {
    fn eq(&self, other: &Self) -> bool {
        self.path == other.path
    }
}

impl<T: FieldBasedMerkleTreeParameters> FieldBasedMerkleTreePath for FieldBasedMHTPath<T> {
    type H = T::H;

    /// A Merkle Path for a leaf of a Merkle Tree with arity >= 2 will be made up of couples of nodes
    /// and an integer. The nodes are all the siblings of the leaf ( in number MERKLE_ARITY - 1 )
    /// and the integer is the position of the leaf, among its siblings, in the input of the hash
    /// function (valid values are from 0 to MERKLE_ARITY - 1).
    type Path = Vec<(Vec<<Self::H as FieldBasedHash>::Data>, usize)>;

    type Parameters = T;

    fn new(path: Self::Path) -> Self {
        Self { path }
    }

    /// NOTE: Check path semantic validity before calling this function.
    fn compute_root(
        &self,
        leaf: &<Self::H as FieldBasedHash>::Data,
    ) -> <Self::H as FieldBasedHash>::Data {
        // Rate may also be smaller than the arity actually, but this assertion
        // is reasonable and simplify the design. Should be also enforced by the
        // MerkleTree that creates this instance, but let's do it again.
        assert_eq!(
            <<Self::H as FieldBasedHash>::Parameters as FieldBasedHashParameters>::R,
            T::MERKLE_ARITY
        );

        let mut digest = <Self::H as FieldBasedHash>::init_constant_length(T::MERKLE_ARITY, None);
        let mut prev_node = leaf.clone();
        for (sibling_nodes, position) in self.path.iter() {
            // Update the digest respecting the position of each sibling
            let mut sibling_idx = 0;
            for i in 0..T::MERKLE_ARITY {
                if i == *position {
                    digest.update(prev_node.clone());
                } else {
                    digest.update(sibling_nodes[sibling_idx]);
                    sibling_idx += 1;
                }
            }

            // Compute the parent node
            prev_node = digest.finalize().unwrap();
            digest.reset(None);
        }

        prev_node
    }

    fn get_raw_path(&self) -> &Self::Path {
        &self.path
    }

    fn get_length(&self) -> usize {
        self.path.len()
    }

    #[inline]
    fn is_leftmost(&self) -> bool {
        for &(_, direction) in &self.path {
            if direction != 0 {
                return false;
            }
        }
        return true;
    }

    #[inline]
    fn is_rightmost(&self) -> bool {
        for &(_, direction) in &self.path {
            if direction != (T::MERKLE_ARITY - 1) {
                return false;
            }
        }
        return true;
    }

    #[inline]
    fn are_right_leaves_empty(&self) -> bool {
        assert!(check_precomputed_parameters::<T>(self.path.len()));

        let mut height = 0usize;
        for &(ref siblings, direction) in &self.path {
            // If the node on the path is not in the rightmost position
            if direction != T::MERKLE_ARITY - 1 {
                // Save the empty node for this height
                let empty_node = T::ZERO_NODE_CST.unwrap().nodes[height].clone();

                // If its following siblings are not the empty nodes, then the node
                // cannot be the non empty rightmost at this height and for the
                // whole tree
                for i in direction..T::MERKLE_ARITY - 1 {
                    if siblings[i] != empty_node {
                        return false;
                    }
                }
            }
            height += 1;
        }
        return true;
    }

    #[inline]
    fn leaf_index(&self) -> usize {
        let mut leaf_index = 0;
        self.path
            .iter()
            .enumerate()
            .for_each(|(i, (_, pos))| leaf_index += T::MERKLE_ARITY.pow(i as u32) * pos);

        leaf_index
    }
}

impl<T: FieldBasedMerkleTreeParameters> ToBytes for FieldBasedMHTPath<T> {
    fn write<W: Write>(&self, mut writer: W) -> IoResult<()> {
        (self.path.len() as u8).write(&mut writer)?;
        for &(ref siblings, position) in self.path.as_slice() {
            siblings.write(&mut writer)?;
            (position as u8).write(&mut writer)?;
        }
        Ok(())
    }
}

impl<T: FieldBasedMerkleTreeParameters> FromBytes for FieldBasedMHTPath<T> {
    fn read<R: Read>(mut reader: R) -> IoResult<Self> {
        let siblings_len = (T::MERKLE_ARITY - 1) as usize;
        let length = u8::read(&mut reader)? as usize;
        let mut path = Vec::with_capacity(length);
        for _ in 0..length {
            let mut siblings = Vec::with_capacity(siblings_len);
            for _ in 0..siblings_len {
                let sibling = <T::H as FieldBasedHash>::Data::read(&mut reader)?;
                siblings.push(sibling);
            }
            let position = u8::read(&mut reader)?;
            path.push((siblings, position as usize));
        }
        Ok(Self { path })
    }
}

/// A wrapper around a Merkle Path for a FieldBasedMerkleTree of arity 2. Merkle Trees of arity
/// 2 are the most common and it's worth to explicitly create a separate struct
#[derive(Derivative)]
#[derivative(
    Clone(bound = ""),
    Debug(bound = ""),
    Default(bound = ""),
    Eq(bound = "")
)]
#[derive(Serialize, Deserialize, CanonicalSerialize, CanonicalDeserialize)]
pub struct FieldBasedBinaryMHTPath<T: FieldBasedMerkleTreeParameters> {
    path: Vec<(<T::H as FieldBasedHash>::Data, bool)>,
}

impl<T: FieldBasedMerkleTreeParameters> SemanticallyValid for FieldBasedBinaryMHTPath<T> {
    fn is_valid(&self) -> bool {
        for (fe, _) in self.path.iter() {
            if !fe.is_valid() {
                return false;
            }
        }
        true
    }
}

impl<T: FieldBasedMerkleTreeParameters> PartialEq for FieldBasedBinaryMHTPath<T> {
    fn eq(&self, other: &Self) -> bool {
        self.path == other.path
    }
}

impl<T: FieldBasedMerkleTreeParameters> FieldBasedMerkleTreePath for FieldBasedBinaryMHTPath<T> {
    type H = T::H;
    type Path = Vec<(<T::H as FieldBasedHash>::Data, bool)>;
    type Parameters = T;

    fn new(path: Self::Path) -> Self {
        Self { path }
    }

    fn compute_root(
        &self,
        leaf: &<Self::H as FieldBasedHash>::Data,
    ) -> <Self::H as FieldBasedHash>::Data {
        // Rate may also be smaller than the arity actually, but this assertion
        // is reasonable and simplify the design. Should be also enforced by the
        // MerkleTree that creates this instance, but let's do it again.
        assert_eq!(
            <<Self::H as FieldBasedHash>::Parameters as FieldBasedHashParameters>::R,
            T::MERKLE_ARITY
        );
        let mut digest = <Self::H as FieldBasedHash>::init_constant_length(2, None);
        let mut prev_node = leaf.clone();
        for (sibling, direction) in self.path.iter() {
            // Choose left and right hash according to direction
            let (left, right) = if !direction {
                (prev_node, sibling.clone())
            } else {
                (sibling.clone(), prev_node)
            };

            // Compute the parent node
            prev_node = digest.update(left).update(right).finalize().unwrap();

            digest.reset(None);
        }

        prev_node
    }

    fn get_raw_path(&self) -> &Self::Path {
        &self.path
    }

    fn get_length(&self) -> usize {
        self.path.len()
    }

    #[inline]
    fn is_leftmost(&self) -> bool {
        for &(_, direction) in &self.path {
            if direction {
                return false;
            }
        }
        return true;
    }

    #[inline]
    fn is_rightmost(&self) -> bool {
        for &(_, direction) in &self.path {
            if !direction {
                return false;
            }
        }
        return true;
    }

    #[inline]
    fn are_right_leaves_empty(&self) -> bool {
        assert!(check_precomputed_parameters::<T>(self.path.len()));

        let mut height = 0usize;
        for &(sibling, direction) in &self.path {
            // If the node on the path is not in the rightmost position
            if !direction {
                // If its following sibling is not the empty node, then the node
                // cannot be the non empty rightmost at this height and for the
                // whole tree
                if sibling != T::ZERO_NODE_CST.unwrap().nodes[height] {
                    return false;
                }
            }
            height += 1;
        }
        return true;
    }

    #[inline]
    fn leaf_index(&self) -> usize {
        let mut leaf_index = 0;
        self.path.iter().enumerate().for_each(|(i, (_, pos))| {
            if *pos {
                leaf_index += 1 << i
            }
        });

        leaf_index as usize
    }
}

/// Serialization utilities for Merkle Path

impl<T: FieldBasedMerkleTreeParameters> ToBytes for FieldBasedBinaryMHTPath<T> {
    fn write<W: Write>(&self, mut writer: W) -> IoResult<()> {
        (self.path.len() as u8).write(&mut writer)?;
        for &(node, direction) in self.path.as_slice() {
            node.write(&mut writer)?;
            direction.write(&mut writer)?;
        }
        Ok(())
    }
}

impl<T: FieldBasedMerkleTreeParameters> FromBytes for FieldBasedBinaryMHTPath<T> {
    fn read<R: Read>(mut reader: R) -> IoResult<Self> {
        let length = u8::read(&mut reader)? as usize;
        let mut path = Vec::with_capacity(length);
        for _ in 0..length {
            let node = <T::H as FieldBasedHash>::Data::read(&mut reader)?;
            let direction = bool::read(&mut reader)?;
            path.push((node, direction));
        }
        Ok(Self { path })
    }
}

/// Conversion utilities for FieldBasedMHTPath and FieldBasedBinaryMHTPath
impl<T: FieldBasedMerkleTreeParameters> From<FieldBasedBinaryMHTPath<T>> for FieldBasedMHTPath<T> {
    fn from(other: FieldBasedBinaryMHTPath<T>) -> Self {
        let mut converted = Vec::with_capacity(other.path.len());
        for &(node, direction) in &other.path {
            converted.push((vec![node], if !direction { 0 } else { 1 }));
        }
        FieldBasedMHTPath::<T>::new(converted)
    }
}

impl<T: FieldBasedMerkleTreeParameters> TryFrom<FieldBasedMHTPath<T>>
    for FieldBasedBinaryMHTPath<T>
{
    type Error = Error;

    fn try_from(other: FieldBasedMHTPath<T>) -> Result<Self, Self::Error> {
        let mut converted = Vec::with_capacity(other.path.len());
        for (nodes, position) in other.path {
            if nodes.len() != 1 {
                Err(format!("There must be only 1 node for each element in the path to be able to perform conversion to a binary path"))?
            }
            if position != 0 && position != 1 {
                Err(format!("Position must be only 0 or 1 for each element in the path to be able to perform conversion to a binary path"))?
            }

            converted.push((nodes[0], if position == 0 { false } else { true }));
        }
        Ok(FieldBasedBinaryMHTPath::<T>::new(converted))
    }
}

impl<T: FieldBasedMerkleTreeParameters> PartialEq<FieldBasedMHTPath<T>>
    for FieldBasedBinaryMHTPath<T>
{
    fn eq(&self, other: &FieldBasedMHTPath<T>) -> bool {
        self == other
    }
}
