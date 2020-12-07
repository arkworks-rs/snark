use crate::{
    crh::*, field_based_mht::*, Error
};
use std::{
    clone::Clone, fmt::Debug, io::{Write, Result as IoResult, Read}
};

/// An implementation of the FieldBasedMerkleTreePath trait, for a given FieldBasedHash and
/// FieldBasedMerkleTree with arbitrary arity.
/// TODO: Test for arity > 2
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct FieldBasedMHTPath<T: FieldBasedMerkleTreeParameters>{
    path: Vec<(Vec<<T::H as FieldBasedHash>::Data>, usize)>,
}

impl<T: FieldBasedMerkleTreeParameters> FieldBasedMHTPath<T> {

    /// Returns true if `self` is a Merkle Path for the left most leaf of a Merkle Tree,
    /// false, otherwise.
    #[inline]
    pub fn is_leftmost(&self) -> bool {
        for &(_, direction) in &self.path {
            if direction != 0 {
                return false;
            }
        }
        return true;
    }

    /// Returns true if `self` is a Merkle Path for the right most leaf of a Merkle Tree,
    /// false, otherwise.
    #[inline]
    pub fn is_rightmost(&self) -> bool {
        for &(_, direction) in &self.path {
            if direction != (T::MERKLE_ARITY - 1) {
                return false;
            }
        }
        return true;
    }

    /// Returns true if `self` is a Merkle Path for the righmost non-empty leaf of the Merkle Tree
    /// (e.g. the leaf which is not physically in the rightmost position of the tree, but it's
    /// followed by all empty leaves).
    /// Assumptions:
    /// 1) Append-only Merkle Tree;
    /// 2) T::EMPTY_HASH_CST is specified;
    /// 3) Not to be called on Merkle Path corresponding to an empty leaf.
    #[inline]
    pub fn is_non_empty_rightmost(&self) -> bool {
        assert!(check_precomputed_parameters::<T>(self.path.len()));

        let mut height = 0usize;
        for &(ref siblings, direction) in &self.path {

            // If the node on the path is not in the rightmost position
            if direction != T::MERKLE_ARITY - 1 {

                // If its following sibling is not the empty node, then the node
                // cannot be the non empty rightmost at this height and for the
                // whole tree
                if siblings[direction] != T::EMPTY_HASH_CST.unwrap().nodes[height] {
                    return false;
                }
            }
            height += 1;
        }
        return true;
    }

    /// Returns the index of the leaf, corresponding to the `self` Merkle Path, in the
    /// corresponding Merkle Tree.
    pub fn leaf_index(&self) -> usize {
        let mut leaf_index = 0;
        self.path
            .iter()
            .enumerate()
            .for_each(|(i, (_, pos))| leaf_index += T::MERKLE_ARITY.pow(i as u32) * pos);

        leaf_index
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

    fn verify_without_length_check(
        &self,
        leaf: &<Self::H as FieldBasedHash>::Data,
        root: &<Self::H as FieldBasedHash>::Data
    ) -> Result<bool, Error> {

        // Rate may also be smaller than the arity actually, but this assertion
        // is reasonable and simplify the design. Should be also enforced by the
        // MerkleTree that creates this instance, but let's do it again.
        assert_eq!(<<Self::H as FieldBasedHash>::Parameters as FieldBasedHashParameters>::R, T::MERKLE_ARITY);

        let mut digest = <Self::H as FieldBasedHash>::init(None);
        let mut prev_node = leaf.clone();
        for (sibling_nodes, position) in self.path.as_slice() {
            assert_eq!(sibling_nodes.len(), T::MERKLE_ARITY - 1);

            // Get the position of the node among its siblings
            let prev_node_position = *position % T::MERKLE_ARITY;

            // Update the digest respecting the position of each sibling
            for i in 0..T::MERKLE_ARITY {
                if i == prev_node_position {
                    digest.update(prev_node.clone());
                } else {
                    let index = i % (T::MERKLE_ARITY - 1); // Be sure to not overflow the siblings vector
                    digest.update(sibling_nodes[index]);
                }
            }

            // Compute the parent node
            prev_node = digest.finalize();
            digest.reset(None);
        }

        // Check final computed node is equal to the root
        if prev_node == root.clone() {
            Ok(true)
        } else {
            Ok(false)
        }
    }

    fn get_raw_path(&self) -> Self::Path {
        self.path.clone()
    }

    fn get_length(&self) -> usize {
        self.path.len()
    }
}

impl<T: FieldBasedMerkleTreeParameters> ToBytes for FieldBasedMHTPath<T> {
    fn write<W: Write>(&self, mut writer: W) -> IoResult<()> {
        (self.path.len() as u8).write(&mut writer)?;
        for &(ref siblings, position) in self.path.as_slice() {
            (siblings.len() as u8).write(&mut writer)?;
            siblings.write(&mut writer)?;
            (position as u8).write(&mut writer)?;
        }
        Ok(())
    }
}

impl<T: FieldBasedMerkleTreeParameters> FromBytes for FieldBasedMHTPath<T> {
    fn read<R: Read>(mut reader: R) -> IoResult<Self> {
        let length = u8::read(&mut reader)? as usize;
        let mut path = Vec::with_capacity(length);
        for _ in 0..length {
            let siblings_len = u8::read(&mut reader)?;
            let mut siblings = Vec::with_capacity(siblings_len as usize);
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
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct FieldBasedBinaryMHTPath<T: FieldBasedMerkleTreeParameters>{
    path: Vec<(<T::H as FieldBasedHash>::Data, bool)>,
}

impl<T: FieldBasedMerkleTreeParameters> FieldBasedBinaryMHTPath<T> {

    /// Returns true if `self` is a Merkle Path for the left most leaf of a Merkle Tree,
    /// false, otherwise.
    #[inline]
    pub fn is_leftmost(&self) -> bool {
        for &(_, direction) in &self.path {
            if direction {
                return false;
            }
        }
        return true;
    }

    /// Returns true if `self` is a Merkle Path for the right most leaf of a Merkle Tree,
    /// false, otherwise.
    #[inline]
    pub fn is_rightmost(&self) -> bool {
        for &(_, direction) in &self.path {
            if !direction {
                return false;
            }
        }
        return true;
    }

    /// Returns true if `self` is a Merkle Path for the righmost non-empty leaf of the Merkle Tree
    /// (e.g. the leaf which is not physically in the rightmost position of the tree, but it's
    /// followed by all empty leaves).
    /// Assumptions:
    /// 1) Append-only Merkle Tree;
    /// 2) T::EMPTY_HASH_CST is specified;
    /// 3) Not to be called on Merkle Path corresponding to an empty leaf.
    #[inline]
    pub fn is_non_empty_rightmost(&self) -> bool {
        assert!(check_precomputed_parameters::<T>(self.path.len()));

        let mut height = 0usize;
        for &(sibling, direction) in &self.path {

            // If the node on the path is not in the rightmost position
            if !direction {

                // If its following sibling is not the empty node, then the node
                // cannot be the non empty rightmost at this height and for the
                // whole tree
                if sibling != T::EMPTY_HASH_CST.unwrap().nodes[height] {
                    return false;
                }
            }
            height += 1;
        }
        return true;
    }

    /// Returns the index of the leaf, corresponding to the `self` Merkle Path, in the
    /// corresponding Merkle Tree.
    pub fn leaf_index(&self) -> usize {
        let mut leaf_index = 0;
        self.path
            .iter()
            .enumerate()
            .for_each(|(i, (_, pos))| {
                if *pos { leaf_index += 1 << i }
            });

        leaf_index as usize
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

    fn verify_without_length_check(
        &self,
        leaf: &<Self::H as FieldBasedHash>::Data,
        root: &<Self::H as FieldBasedHash>::Data
    ) -> Result<bool, Error> {

        // Rate may also be smaller than the arity actually, but this assertion
        // is reasonable and simplify the design. Should be also enforced by the
        // MerkleTree that creates this instance, but let's do it again.
        assert_eq!(<<Self::H as FieldBasedHash>::Parameters as FieldBasedHashParameters>::R, T::MERKLE_ARITY);
        let mut digest = <Self::H as FieldBasedHash>::init(None);
        let mut prev_node = *leaf;
        for &(sibling, direction) in self.path.as_slice() {

            // Choose left and right hash according to direction
            let (left, right) = if !direction {
                (prev_node, sibling)
            } else {
                (sibling, prev_node)
            };

            // Compute the parent node
            prev_node = digest
                .update(left)
                .update(right)
                .finalize();

            digest.reset(None);
        }

        // Check final computed node is equal to the root
        if prev_node == root.clone() {
            Ok(true)
        } else {
            Ok(false)
        }
    }

    fn get_raw_path(&self) -> Self::Path {
        self.path.clone()
    }

    fn get_length(&self) -> usize {
        self.path.len()
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
            converted.push((vec![node], if !direction {0} else {1}));
        }
        FieldBasedMHTPath::<T>::new(converted)
    }
}

impl<T: FieldBasedMerkleTreeParameters> From<FieldBasedMHTPath<T>> for FieldBasedBinaryMHTPath<T> {
    fn from(other: FieldBasedMHTPath<T>) -> Self {
        let mut converted = Vec::with_capacity(other.path.len());
        for (nodes, position) in other.path {
            assert!(nodes.len() == 1);
            assert!(position == 0 || position == 1);

            converted.push((nodes[0], if position == 0 {false} else {true}));
        }
        FieldBasedBinaryMHTPath::<T>::new(converted)
    }
}

impl<T: FieldBasedMerkleTreeParameters> PartialEq<FieldBasedMHTPath<T>> for FieldBasedBinaryMHTPath<T> {
    fn eq(&self, other: &FieldBasedMHTPath<T>) -> bool {
        self == other
    }
}