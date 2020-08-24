// Implements AddAssign on Self by deferring to an implementation on &Self
#[macro_export]
macro_rules! impl_ops {
    (
        $type: ty,
        $native: ty,
        $trait: ident,
        $fn: ident,
        $assign_trait: ident,
        $assign_fn: ident,
        $impl: expr,
        $constant_impl: expr,
        $($args:tt)*
    ) => {
        impl_bounded_ops!($type, $native, $trait, $fn, $assign_trait, $assign_fn, $impl, $constant_impl, ($($args)+), );
    };
}

macro_rules! impl_bounded_ops {
    (
        $type: ty,
        $native: ty,
        $trait: ident,
        $fn: ident,
        $assign_trait: ident,
        $assign_fn: ident,
        $impl: expr,
        $constant_impl: expr,
        ($($params:tt)+),
        $($bounds:tt)*
    ) => {
        impl<'a, $($params)+> core::ops::$trait<&'a $type> for &'a $type
        where
            $($bounds)*
        {
            type Output = $type;

            fn $fn(self, other: Self) -> Self::Output {
                $impl(self, other)
            }
        }

        impl<'a, $($params)+> core::ops::$trait<$type> for &'a $type
        where
            $($bounds)*
        {
            type Output = $type;

            fn $fn(self, other: $type) -> Self::Output {
                core::ops::$trait::$fn(self, &other)
            }
        }

        impl<'a, $($params)+> core::ops::$trait<&'a $type> for $type
        where
            $($bounds)*
        {
            type Output = $type;

            fn $fn(self, other: &'a $type) -> Self::Output {
                core::ops::$trait::$fn(&self, other)
            }
        }

        impl<$($params)+> core::ops::$trait<$type> for $type
        where

            $($bounds)*
        {
            type Output = $type;

            fn $fn(self, other: $type) -> Self::Output {
                core::ops::$trait::$fn(&self, &other)
            }
        }

        impl<$($params)+> core::ops::$assign_trait<$type> for $type
        where

            $($bounds)*
        {
            fn $assign_fn(&mut self, other: $type) {
                let result = core::ops::$trait::$fn(&*self, &other);
                *self = result
            }
        }

        impl<'a, $($params)+> core::ops::$assign_trait<&'a $type> for $type
        where

            $($bounds)*
        {
            fn $assign_fn(&mut self, other: &'a $type) {
                let result = core::ops::$trait::$fn(&*self, other);
                *self = result
            }
        }

        impl<'a, $($params)+> core::ops::$trait<$native> for &'a $type
        where

            $($bounds)*
        {
            type Output = $type;

            fn $fn(self, other: $native) -> Self::Output {
                $constant_impl(self, other)
            }
        }

        impl<$($params)+> core::ops::$trait<$native> for $type
        where

            $($bounds)*
        {
            type Output = $type;

            fn $fn(self, other: $native) -> Self::Output {
                core::ops::$trait::$fn(&self, other)
            }
        }

        impl<$($params)+> core::ops::$assign_trait<$native> for $type
        where

            $($bounds)*
        {

            fn $assign_fn(&mut self, other: $native) {
                let result = core::ops::$trait::$fn(&*self, other);
                *self = result
            }
        }
    }
}
