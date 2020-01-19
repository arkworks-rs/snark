extern crate proc_macro;

use proc_macro::TokenStream;
use quote::quote;
use syn::{parse_macro_input, DeriveInput};

// This assumes you already have `Add` defined on `&'a Self` operand
#[proc_macro_derive(AddFromRef)]
pub fn derive_add_from_ref(input: TokenStream) -> TokenStream {
    let item = parse_macro_input!(input as DeriveInput);

    let name = &item.ident;
    let (impl_generics, ty_generics, where_clause) = item.generics.split_for_impl();

    quote!(
        impl #impl_generics ::std::ops::Add<Self> for #name #ty_generics
        #where_clause
        {
            type Output = Self;

            #[inline]
            fn add(self, other: Self) -> Self {
                <Self as Add<&Self>>::add(self, &other)
            }
        }
    )
    .into()
}

// This assumes you already have `AddAssign` defined on a `&'a Self` operand
#[proc_macro_derive(AddAssignFromRef)]
pub fn derive_addassign_from_ref(input: TokenStream) -> TokenStream {
    let item = parse_macro_input!(input as DeriveInput);

    let name = &item.ident;
    let (impl_generics, ty_generics, where_clause) = item.generics.split_for_impl();

    quote!(
        impl #impl_generics ::std::ops::AddAssign<Self> for #name #ty_generics
        #where_clause
        {
            fn add_assign(&mut self, other: Self) {
                <Self as ::std::ops::AddAssign<&Self>>::add_assign(self, &other)
            }
        }
    )
    .into()
}

// This assumes you already have `Sub` defined on `&'a Self` operand
#[proc_macro_derive(SubFromRef)]
pub fn derive_sub_from_ref(input: TokenStream) -> TokenStream {
    let item = parse_macro_input!(input as DeriveInput);

    let name = &item.ident;
    let (impl_generics, ty_generics, where_clause) = item.generics.split_for_impl();

    quote!(
        impl #impl_generics ::std::ops::Sub<Self> for #name #ty_generics
        #where_clause
        {
            type Output = Self;

            #[inline]
            fn sub(self, other: Self) -> Self {
                <Self as Sub<&Self>>::sub(self, &other)
            }
        }
    )
    .into()
}

// This assumes you already have `SubAssign` defined on a `&'a Self` operand
#[proc_macro_derive(SubAssignFromRef)]
pub fn derive_subassign_from_ref(input: TokenStream) -> TokenStream {
    let item = parse_macro_input!(input as DeriveInput);

    let name = &item.ident;
    let (impl_generics, ty_generics, where_clause) = item.generics.split_for_impl();

    quote!(
        impl #impl_generics ::std::ops::SubAssign<Self> for #name #ty_generics
        #where_clause
        {
            fn sub_assign(&mut self, other: Self) {
                <Self as ::std::ops::SubAssign<&Self>>::sub_assign(self, &other)
            }
        }
    )
    .into()
}

// This assumes you already have `Mul` defined on `&'a Self` operand
#[proc_macro_derive(MulFromRef)]
pub fn derive_mul_from_ref(input: TokenStream) -> TokenStream {
    let item = parse_macro_input!(input as DeriveInput);

    let name = &item.ident;
    let (impl_generics, ty_generics, where_clause) = item.generics.split_for_impl();

    quote!(
        impl #impl_generics ::std::ops::Mul<Self> for #name #ty_generics
        #where_clause
        {
            type Output = Self;

            #[inline]
            fn mul(self, other: Self) -> Self {
                <Self as Mul<&Self>>::mul(self, &other)
            }
        }
    )
    .into()
}

// This assumes you already have `MulAssign` defined on a `&'a Self` operand
#[proc_macro_derive(MulAssignFromRef)]
pub fn derive_mulassign_from_ref(input: TokenStream) -> TokenStream {
    let item = parse_macro_input!(input as DeriveInput);

    let name = &item.ident;
    let (impl_generics, ty_generics, where_clause) = item.generics.split_for_impl();

    quote!(
        impl #impl_generics ::std::ops::MulAssign<Self> for #name #ty_generics
        #where_clause
        {
            fn mul_assign(&mut self, other: Self) {
                <Self as ::std::ops::MulAssign<&Self>>::mul_assign(self, &other)
            }
        }
    )
    .into()
}

// This assumes you already have `Div` defined on `&'a Self` operand
#[proc_macro_derive(DivFromRef)]
pub fn derive_div_from_ref(input: TokenStream) -> TokenStream {
    let item = parse_macro_input!(input as DeriveInput);

    let name = &item.ident;
    let (impl_generics, ty_generics, where_clause) = item.generics.split_for_impl();

    quote!(
        impl #impl_generics ::std::ops::Div<Self> for #name #ty_generics
        #where_clause
        {
            type Output = Self;

            #[inline]
            fn div(self, other: Self) -> Self {
                <Self as Div<&Self>>::div(self, &other)
            }
        }
    )
    .into()
}

// This assumes you already have `DivAssign` defined on a `&'a Self` operand
#[proc_macro_derive(DivAssignFromRef)]
pub fn derive_divassign_from_ref(input: TokenStream) -> TokenStream {
    let item = parse_macro_input!(input as DeriveInput);

    let name = &item.ident;
    let (impl_generics, ty_generics, where_clause) = item.generics.split_for_impl();

    quote!(
        impl #impl_generics ::std::ops::DivAssign<Self> for #name #ty_generics
        #where_clause
        {
            fn div_assign(&mut self, other: Self) {
                <Self as ::std::ops::DivAssign<&Self>>::div_assign(self, &other)
            }
        }
    )
    .into()
}
