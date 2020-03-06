extern crate proc_macro;

use proc_macro2::TokenStream;
use syn::{parse_macro_input, Data, DeriveInput, Index, Type};

use quote::{quote, ToTokens};

#[proc_macro_derive(CanonicalSerialize)]
pub fn derive_canonical_serialize(input: proc_macro::TokenStream) -> proc_macro::TokenStream {
    let ast = parse_macro_input!(input as DeriveInput);
    proc_macro::TokenStream::from(impl_canonical_serialize(&ast))
}

fn impl_serialize_field(
    serialize_body: &mut Vec<TokenStream>,
    serialized_size_body: &mut Vec<TokenStream>,
    idents: &mut Vec<Box<dyn ToTokens>>,
    ty: &Type,
) {
    // Check if type is a tuple.
    match ty {
        Type::Tuple(tuple) => {
            for (i, elem_ty) in tuple.elems.iter().enumerate() {
                let index = Index::from(i);
                idents.push(Box::new(index));
                impl_serialize_field(serialize_body, serialized_size_body, idents, elem_ty);
                idents.pop();
            }
        },
        _ => {
            serialize_body
                .push(quote! { CanonicalSerialize::serialize(&self.#(#idents).*, writer)?; });
            serialized_size_body
                .push(quote! { size += CanonicalSerialize::serialized_size(&self.#(#idents).*); });
        },
    }
}

fn impl_canonical_serialize(ast: &syn::DeriveInput) -> TokenStream {
    let name = &ast.ident;

    let (impl_generics, ty_generics, where_clause) = ast.generics.split_for_impl();

    let mut serialize_body = Vec::<TokenStream>::new();
    let mut serialized_size_body = Vec::<TokenStream>::new();

    match ast.data {
        Data::Struct(ref data_struct) => {
            for (i, field) in data_struct.fields.iter().enumerate() {
                let mut idents = Vec::<Box<dyn ToTokens>>::new();
                match field.ident {
                    None => {
                        let index = Index::from(i);
                        idents.push(Box::new(index));
                    },
                    Some(ref ident) => {
                        idents.push(Box::new(ident.clone()));
                    },
                }

                impl_serialize_field(
                    &mut serialize_body,
                    &mut serialized_size_body,
                    &mut idents,
                    &field.ty,
                );
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
            fn serialize<W: Write>(&self, writer: &mut W) -> Result<(), SerializationError> {
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

fn impl_deserialize_field(ty: &Type) -> TokenStream {
    // Check if type is a tuple.
    match ty {
        Type::Tuple(tuple) => {
            let mut fields = Vec::new();
            for elem_ty in tuple.elems.iter() {
                fields.push(impl_deserialize_field(elem_ty));
            }
            quote! { (#(#fields)*), }
        },
        _ => {
            quote! { CanonicalDeserialize::deserialize(reader)?, }
        },
    }
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
                        field_cases.push(impl_deserialize_field(&field.ty))
                    },
                    // struct field without len_type
                    Some(ident) => {
                        let field = impl_deserialize_field(&field.ty);
                        field_cases.push(quote! { #ident: #field })
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
            fn deserialize<R: Read>(reader: &mut R) -> Result<Self, SerializationError> {
                #deserialize_body
            }
        }
    };
    gen
}
