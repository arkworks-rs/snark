extern crate proc_macro;

use proc_macro::TokenStream;
use quote::quote;
use syn::{parse_macro_input, DeriveInput};

fn get_type_from_attrs(attrs: &[syn::Attribute], attr_name: &str) -> syn::Result<syn::LitStr> {
    attrs
        .iter()
        .find(|attr| attr.path.is_ident(attr_name))
        .map_or_else(
            || {
                Err(syn::Error::new(
                    proc_macro2::Span::call_site(),
                    format!("Could not find attribute {}", attr_name),
                ))
            },
            |attr| match attr.parse_meta()? {
                syn::Meta::NameValue(meta) => {
                    if let syn::Lit::Str(lit) = &meta.lit {
                        Ok(lit.clone())
                    } else {
                        Err(syn::Error::new_spanned(
                            meta,
                            &format!("Could not parse {} attribute", attr_name)[..],
                        ))
                    }
                },
                bad => Err(syn::Error::new_spanned(
                    bad,
                    &format!("Could not parse {} attribute", attr_name)[..],
                )),
            },
        )
}

fn get_param_from_generics(
    item: &DeriveInput,
) -> Result<Vec<proc_macro2::TokenStream>, syn::Error> {
    item.generics
        .params
        .iter()
        .map(|generic_param| match generic_param {
            syn::GenericParam::Type(type_param) => {
                let ident = type_param.ident.clone();
                Ok(quote!(#ident))
            },
            bad => Err(syn::Error::new_spanned(
                bad,
                &format!("Could not parse generics of {}", item.ident)[..],
            )),
        })
        .collect()
}

// This assumes you already have `Add` defined on `&'a Self` operand
#[proc_macro_derive(AddFromRef, attributes(ArithmeticBound))]
pub fn derive_add_from_ref(input: TokenStream) -> TokenStream {
    let item = parse_macro_input!(input as DeriveInput);
    let name = &item.ident;

    let (_, ty_generics, where_clause) = item.generics.split_for_impl();
    let bounded_generics: Result<Vec<_>, _> = get_param_from_generics(&item);
    let first_param = {
        let vec = bounded_generics.expect("Could not find generic parameter");

        vec.get(0)
            .expect("Struct should have a generic parameter")
            .clone()
    };

    let arithmetic_bound = get_type_from_attrs(&item.attrs, "ArithmeticBound").unwrap();
    let bound: syn::Type = arithmetic_bound.parse().unwrap();

    quote!(
        impl <#first_param: #bound> ::std::ops::Add<Self> for #name #ty_generics
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
#[proc_macro_derive(AddAssignFromRef, attributes(ArithmeticBound))]
pub fn derive_addassign_from_ref(input: TokenStream) -> TokenStream {
    let item = parse_macro_input!(input as DeriveInput);
    let name = &item.ident;

    let (_, ty_generics, where_clause) = item.generics.split_for_impl();
    let bounded_generics: Result<Vec<_>, _> = get_param_from_generics(&item);
    let first_param = {
        let vec = bounded_generics.expect("Could not find generic parameter");

        vec.get(0)
            .expect("Struct should have a generic parameter")
            .clone()
    };

    let arithmetic_bound = get_type_from_attrs(&item.attrs, "ArithmeticBound").unwrap();
    let bound: syn::Type = arithmetic_bound.parse().unwrap();

    quote!(
        impl <#first_param: #bound> ::std::ops::AddAssign<Self> for #name #ty_generics
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
#[proc_macro_derive(SubFromRef, attributes(ArithmeticBound))]
pub fn derive_sub_from_ref(input: TokenStream) -> TokenStream {
    let item = parse_macro_input!(input as DeriveInput);
    let name = &item.ident;

    let (_, ty_generics, where_clause) = item.generics.split_for_impl();
    let bounded_generics: Result<Vec<_>, _> = get_param_from_generics(&item);
    let first_param = {
        let vec = bounded_generics.expect("Could not find generic parameter");

        vec.get(0)
            .expect("Struct should have a generic parameter")
            .clone()
    };

    let arithmetic_bound = get_type_from_attrs(&item.attrs, "ArithmeticBound").unwrap();
    let bound: syn::Type = arithmetic_bound.parse().unwrap();

    quote!(
        impl <#first_param: #bound> ::std::ops::Sub<Self> for #name #ty_generics
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
#[proc_macro_derive(SubAssignFromRef, attributes(ArithmeticBound))]
pub fn derive_subassign_from_ref(input: TokenStream) -> TokenStream {
    let item = parse_macro_input!(input as DeriveInput);
    let name = &item.ident;

    let (_, ty_generics, where_clause) = item.generics.split_for_impl();
    let bounded_generics: Result<Vec<_>, _> = get_param_from_generics(&item);
    let first_param = {
        let vec = bounded_generics.expect("Could not find generic parameter");

        vec.get(0)
            .expect("Struct should have a generic parameter")
            .clone()
    };

    let arithmetic_bound = get_type_from_attrs(&item.attrs, "ArithmeticBound").unwrap();
    let bound: syn::Type = arithmetic_bound.parse().unwrap();

    quote!(
        impl <#first_param: #bound> ::std::ops::SubAssign<Self> for #name #ty_generics
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
#[proc_macro_derive(MulFromRef, attributes(ArithmeticBound))]
pub fn derive_mul_from_ref(input: TokenStream) -> TokenStream {
    let item = parse_macro_input!(input as DeriveInput);
    let name = &item.ident;

    let (_, ty_generics, where_clause) = item.generics.split_for_impl();
    let bounded_generics: Result<Vec<_>, _> = get_param_from_generics(&item);
    let first_param = {
        let vec = bounded_generics.expect("Could not find generic parameter");

        vec.get(0)
            .expect("Struct should have a generic parameter")
            .clone()
    };

    let arithmetic_bound = get_type_from_attrs(&item.attrs, "ArithmeticBound").unwrap();
    let bound: syn::Type = arithmetic_bound.parse().unwrap();

    quote!(
        impl <#first_param: #bound> ::std::ops::Mul<Self> for #name #ty_generics
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
#[proc_macro_derive(MulAssignFromRef, attributes(ArithmeticBound))]
pub fn derive_mulassign_from_ref(input: TokenStream) -> TokenStream {
    let item = parse_macro_input!(input as DeriveInput);
    let name = &item.ident;

    let (_, ty_generics, where_clause) = item.generics.split_for_impl();
    let bounded_generics: Result<Vec<_>, _> = get_param_from_generics(&item);
    let first_param = {
        let vec = bounded_generics.expect("Could not find generic parameter");

        vec.get(0)
            .expect("Struct should have a generic parameter")
            .clone()
    };

    let arithmetic_bound = get_type_from_attrs(&item.attrs, "ArithmeticBound").unwrap();
    let bound: syn::Type = arithmetic_bound.parse().unwrap();

    quote!(
        impl <#first_param: #bound> ::std::ops::MulAssign<Self> for #name #ty_generics
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
#[proc_macro_derive(DivFromRef, attributes(ArithmeticBound))]
pub fn derive_div_from_ref(input: TokenStream) -> TokenStream {
    let item = parse_macro_input!(input as DeriveInput);
    let name = &item.ident;

    let (_, ty_generics, where_clause) = item.generics.split_for_impl();
    let bounded_generics: Result<Vec<_>, _> = get_param_from_generics(&item);
    let first_param = {
        let vec = bounded_generics.expect("Could not find generic parameter");

        vec.get(0)
            .expect("Struct should have a generic parameter")
            .clone()
    };

    let arithmetic_bound = get_type_from_attrs(&item.attrs, "ArithmeticBound").unwrap();
    let bound: syn::Type = arithmetic_bound.parse().unwrap();

    quote!(
        impl <#first_param: #bound> ::std::ops::Div<Self> for #name #ty_generics
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
#[proc_macro_derive(DivAssignFromRef, attributes(ArithmeticBound))]
pub fn derive_divassign_from_ref(input: TokenStream) -> TokenStream {
    let item = parse_macro_input!(input as DeriveInput);
    let name = &item.ident;

    let (_, ty_generics, where_clause) = item.generics.split_for_impl();
    let bounded_generics: Result<Vec<_>, _> = get_param_from_generics(&item);
    let first_param = {
        let vec = bounded_generics.expect("Could not find generic parameter");

        vec.get(0)
            .expect("Struct should have a generic parameter")
            .clone()
    };

    let arithmetic_bound = get_type_from_attrs(&item.attrs, "ArithmeticBound").unwrap();
    let bound: syn::Type = arithmetic_bound.parse().unwrap();

    quote!(
        impl <#first_param: #bound> ::std::ops::DivAssign<Self> for #name #ty_generics
        #where_clause
        {
            fn div_assign(&mut self, other: Self) {
                <Self as ::std::ops::DivAssign<&Self>>::div_assign(self, &other)
            }
        }
    )
    .into()
}
