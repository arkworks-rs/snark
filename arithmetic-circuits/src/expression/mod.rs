use crate::arithmetic_circuit::{filter_constants, ArithmeticCircuit, Node};
use ark_ff::PrimeField;
use ark_std::{
    collections::HashMap,
    fmt::Display,
    iter::{Product, Sum},
    ops::{Add, AddAssign, Mul, MulAssign, Neg, Sub, SubAssign},
    rc::Rc,
};
use itertools::Itertools;

#[cfg(any(feature = "examples", test))]
pub mod examples;

#[cfg(test)]
mod tests;

/// Utilities that expose a user-friendly way to construct arithmetic circuits,
/// with syntax along the lines of:
/// let x = Expression::Variable("x");
/// let y = Expression::Variable("y");
/// let output = y.pow(2) - x.pow(3) + 1

/// Syntax summary:
/// - Expression::variable(id) creates a variable with the given ID.
///
/// - Expression::constant(value) creates a constant with the given F value.
///
/// - +, - and * are overloaded to mean addition, subtraction a nd multiplication of expressions
///   Their assigning counterparts +=, -=, *= are also overloaded.
///
/// - Constants in the form of F can be used as operands on the right-hand side only.
///   This is due to the implementation for i32 from the next point.
///   E.g.: F::from(3) * exp, F::ONE * exp, and exp * F::from(3) are all valid
///   However, 3 * exp, -5 * exp, and exp * 3 are not.
///
/// - Constants in the form of i32 (where F: From<i32>) can be used as operands on the left-hand side only.
///   This is due to i32 and PrimeField both being foreign types.
///   E.g. 1 + exp and -5 * exp are both valid, equivalent to F::from(1) + exp and F::from(-5) * exp, respectively.
///   However, exp + 1, exp - 3 and exp * -5 are not.
enum ExpressionInner<F: PrimeField> {
    Variable(String),
    Constant(F),
    Add(Expression<F>, Expression<F>),
    Mul(Expression<F>, Expression<F>),
}

// New type pattern necessary so that we can implement operators such as +,
// which we can't directly do on the foreign type Rc<ExpressionInner<F>>
pub struct Expression<F: PrimeField>(Rc<ExpressionInner<F>>);

impl<F: PrimeField> Expression<F> {
    pub fn constant(value: F) -> Self {
        Expression(Rc::new(ExpressionInner::Constant(value)))
    }

    pub fn variable(label: &str) -> Self {
        Expression(Rc::new(ExpressionInner::Variable(label.to_string())))
    }

    pub fn to_arithmetic_circuit(&self) -> ArithmeticCircuit<F> {
        let mut nodes = HashMap::new();
        self.update_map(&mut nodes);

        let ptr_to_idx = nodes
            .iter()
            .map(|(ptr, (idx, _))| (*ptr, nodes.len() - idx - 1))
            .collect::<HashMap<_, _>>();

        let sorted_nodes = nodes
            .into_iter()
            .sorted_by(|(_, (i, _)), (_, (j, _))| j.cmp(i))
            .map(|(_, (_, node))| node)
            .collect::<Vec<_>>();

        let mut nodes = Vec::new();
        for node in sorted_nodes {
            match node {
                Node::Variable(label) => {
                    nodes.push(Node::Variable(label));
                }
                Node::Constant(value) => {
                    nodes.push(Node::Constant(value));
                }
                Node::Add(a, b) => {
                    nodes.push(Node::Add(ptr_to_idx[&a], ptr_to_idx[&b]));
                }
                Node::Mul(a, b) => {
                    nodes.push(Node::Mul(ptr_to_idx[&a], ptr_to_idx[&b]));
                }
            }
        }

        let (nodes, constants) = filter_constants(&nodes);

        let variables = HashMap::from_iter(nodes.iter().enumerate().filter_map(|(i, node)| {
            if let Node::Variable(label) = node {
                Some((label.clone(), i))
            } else {
                None
            }
        }));

        ArithmeticCircuit {
            nodes,
            constants,
            variables,
            unit_group_bits: None,
        }
    }

    fn pointer(&self) -> usize {
        self.0.as_ref() as *const _ as usize
    }

    fn update_map(&self, nodes: &mut HashMap<usize, (usize, Node<F>)>) {
        if nodes.contains_key(&self.pointer()) {
            return;
        }
        match &*self.0 {
            ExpressionInner::Variable(label) => {
                nodes.insert(self.pointer(), (nodes.len(), Node::Variable(label.clone())));
            }
            ExpressionInner::Constant(value) => {
                nodes.insert(self.pointer(), (nodes.len(), Node::Constant(*value)));
            }
            ExpressionInner::Add(a, b) => {
                nodes.insert(
                    self.pointer(),
                    (nodes.len(), Node::Add(a.pointer(), b.pointer())),
                );
                a.update_map(nodes);
                b.update_map(nodes);
            }
            ExpressionInner::Mul(a, b) => {
                nodes.insert(
                    self.pointer(),
                    (nodes.len(), Node::Mul(a.pointer(), b.pointer())),
                );
                a.update_map(nodes);
                b.update_map(nodes);
            }
        }
    }

    pub fn scalar_product(a: Vec<Expression<F>>, b: Vec<Expression<F>>) -> Expression<F> {
        a.into_iter().zip(b).map(|(a, b)| a * b).sum()
    }

    pub fn sparse_scalar_product(a: &Vec<(F, usize)>, b: &Vec<Expression<F>>) -> Expression<F> {
        a.iter()
            .map(|(a, i)| b[*i].clone() * *a)
            .collect::<Vec<_>>()
            .into_iter()
            .sum()
    }

    pub fn pow(self, rhs: usize) -> Self {
        if rhs == 0 {
            return self;
        }

        let mut bits = (0..usize::BITS).rev().map(|pos| (rhs >> pos) & 1);

        bits.position(|bit| bit == 1);

        let mut current = self.clone();

        for bit in bits {
            current = current.clone() * current;

            if bit == 1 {
                current = current.clone() * self.clone();
            }
        }

        current
    }
}

impl<F: PrimeField> Clone for Expression<F> {
    fn clone(&self) -> Self {
        Expression(Rc::clone(&self.0))
    }
}

impl<F: PrimeField> Neg for Expression<F> {
    type Output = Expression<F>;

    fn neg(self) -> Self::Output {
        Expression::constant(-F::ONE) * self
    }
}

impl<F: PrimeField> Add for Expression<F> {
    type Output = Expression<F>;

    fn add(self, rhs: Expression<F>) -> Self::Output {
        Expression(Rc::new(ExpressionInner::Add(self.clone(), rhs.clone())))
    }
}

impl<F: PrimeField> Mul for Expression<F> {
    type Output = Expression<F>;

    fn mul(self, rhs: Self) -> Self::Output {
        Expression(Rc::new(ExpressionInner::Mul(self.clone(), rhs.clone())))
    }
}

impl<F: PrimeField> Sub for Expression<F> {
    type Output = Expression<F>;

    fn sub(self, rhs: Self) -> Self::Output {
        self + (-rhs)
    }
}

macro_rules! impl_constant_op {
    ($op_trait:ident, $method:ident, $op:tt) => {
        impl<F: PrimeField + From<i32>> $op_trait<Expression<F>> for i32 {
            type Output = Expression<F>;

            fn $method(self, rhs: Expression<F>) -> Self::Output {
                Expression::constant(F::from(self)) $op rhs
            }
        }

        impl<F: PrimeField> $op_trait<F> for Expression<F> {
            type Output = Expression<F>;

            fn $method(self, rhs: F) -> Self::Output {
                self $op Expression::constant(rhs)
            }
        }
    };
}

impl_constant_op!(Add, add, +);
impl_constant_op!(Mul, mul, *);
impl_constant_op!(Sub, sub, -);

macro_rules! impl_op_assign_aux {
    ($op_trait_assign:ident, $method_assign:ident, $op:tt, $self_type:ty) => {
        impl<F: PrimeField> $op_trait_assign<$self_type> for Expression<F> {
            fn $method_assign(&mut self, rhs: $self_type) {
                *self = self.clone() $op rhs;
            }
        }
    };
}

macro_rules! impl_op_assign {
    ($op_trait_assign:ident, $method_assign:ident, $op:tt) => {
        impl_op_assign_aux!($op_trait_assign, $method_assign, $op, F);
        impl_op_assign_aux!($op_trait_assign, $method_assign, $op, Expression<F>);
    };
}

impl_op_assign!(AddAssign, add_assign, +);
impl_op_assign!(MulAssign, mul_assign, *);
impl_op_assign!(SubAssign, sub_assign, -);

impl<F: PrimeField> Sum for Expression<F> {
    fn sum<I: Iterator<Item = Self>>(iter: I) -> Self {
        iter.reduce(|a, b| a + b).unwrap()
    }
}

impl<F: PrimeField> Product for Expression<F> {
    fn product<I: Iterator<Item = Self>>(iter: I) -> Self {
        iter.reduce(|a, b| a * b).unwrap()
    }
}

impl<F: PrimeField> Display for Expression<F> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let hash = self.pointer();
        match &*self.0 {
            ExpressionInner::Variable(label) => write!(f, "Variable({})<{}>", label, hash),
            ExpressionInner::Constant(value) => write!(f, "Constant({:?})<{}>", value, hash),
            ExpressionInner::Add(a, b) => {
                write!(f, "Add({}, {})<{}>", a.pointer(), b.pointer(), hash)
            }
            ExpressionInner::Mul(a, b) => {
                write!(f, "Mul({}, {})<{}>", a.pointer(), b.pointer(), hash)
            }
        }
    }
}
