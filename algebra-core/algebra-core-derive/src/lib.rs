extern crate proc_macro;

use proc_macro2::TokenStream;
use syn::{parse_macro_input, Data, DeriveInput, Index};

use quote::quote;

#[proc_macro_derive(CanonicalSerialize)]
pub fn derive_canonical_serialize(input: proc_macro::TokenStream) -> proc_macro::TokenStream {
    let ast = parse_macro_input!(input as DeriveInput);
    proc_macro::TokenStream::from(impl_canonical_serialize(&ast))
}

fn impl_canonical_serialize(ast: &syn::DeriveInput) -> TokenStream {
    let name = &ast.ident;

    let (impl_generics, ty_generics, where_clause) = ast.generics.split_for_impl();

    let mut serialize_body = Vec::<TokenStream>::new();
    let mut serialized_size_body = Vec::<TokenStream>::new();

    match ast.data {
        Data::Struct(ref data_struct) => {
            for (i, field) in data_struct.fields.iter().enumerate() {
                match field.ident {
                    None => {
                        let index = Index::from(i);
                        serialize_body
                            .push(quote! { CanonicalSerialize::serialize(&self.#index, writer)?; });
                        serialized_size_body.push(
                            quote! { size += CanonicalSerialize::serialized_size(&self.#index); },
                        );
                    },
                    Some(ref ident) => {
                        serialize_body
                            .push(quote! { CanonicalSerialize::serialize(&self.#ident, writer)?; });
                        serialized_size_body.push(
                            quote! { size += CanonicalSerialize::serialized_size(&self.#ident); },
                        );
                    },
                }
            }
        },
        _ => panic!(
            "Serialize can only be derived for structs, {} is not a struct",
            name
        ),
    };

    let gen = quote! {
        impl #impl_generics CanonicalSerialize for #name #ty_generics #where_clause {
            #[allow(unused_mut, unused_variables)]
            fn serialize<W: ::algebra_core::io::Write>(&self, writer: &mut W) -> Result<(), ::algebra_core::SerializationError> {
                #(#serialize_body)*
                Ok(())
            }
            #[allow(unused_mut, unused_variables)]
            fn serialized_size(&self) -> usize {
                let mut size = 0;
                #(#serialized_size_body)*
                size
            }
        }
    };
    gen
}

#[proc_macro_derive(CanonicalDeserialize)]
pub fn derive_canonical_deserialize(input: proc_macro::TokenStream) -> proc_macro::TokenStream {
    let ast = parse_macro_input!(input as DeriveInput);
    proc_macro::TokenStream::from(impl_canonical_deserialize(&ast))
}

fn impl_canonical_deserialize(ast: &syn::DeriveInput) -> TokenStream {
    let name = &ast.ident;

    let (impl_generics, ty_generics, where_clause) = ast.generics.split_for_impl();

    let deserialize_body;

    match ast.data {
        Data::Struct(ref data_struct) => {
            let mut tuple = false;
            let mut field_cases = Vec::<TokenStream>::new();
            for field in data_struct.fields.iter() {
                match &field.ident {
                    None => {
                        tuple = true;
                        field_cases.push(quote! { ::algebra_core::CanonicalDeserialize::deserialize(reader)?, })
                    },
                    // struct field without len_type
                    Some(ident) => {
                        field_cases.push(quote! { #ident: ::algebra_core::CanonicalDeserialize::deserialize(reader)?, })
                    },
                }
            }

            if tuple {
                deserialize_body = quote!({
                    Ok(#name (
                        #(#field_cases)*
                    ))
                });
            } else {
                deserialize_body = quote!({
                    Ok(#name {
                        #(#field_cases)*
                    })
                });
            }
        },
        _ => panic!(
            "Deserialize can only be derived for structs, {} is not a Struct",
            name
        ),
    };

    let gen = quote! {
        impl #impl_generics CanonicalDeserialize for #name #ty_generics #where_clause {
            #[allow(unused_mut,unused_variables)]
            fn deserialize<R: ::algebra_core::io::Read>(reader: &mut R) -> Result<Self, ::algebra_core::SerializationError> {
                #deserialize_body
            }
        }
    };
    gen
}
