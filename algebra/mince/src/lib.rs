#![recursion_limit = "256"]
#![allow(
    clippy::upper_case_acronyms,
    clippy::too_many_arguments,
    clippy::type_complexity,
    clippy::try_err,
    clippy::map_collect_result_unit,
    clippy::not_unsafe_ptr_arg_deref,
    clippy::suspicious_op_assign_impl,
    clippy::suspicious_arithmetic_impl,
    clippy::assertions_on_constants
)]

extern crate proc_macro;

mod intrinsics;
use intrinsics::*;

mod arithmetic;
use arithmetic::*;

use proc_macro::TokenStream;
use quote::quote;

#[proc_macro_attribute]
pub fn assemble(_meta: TokenStream, input: TokenStream) -> TokenStream {
    let ast: syn::ItemFn = syn::parse(input).unwrap();
    let sig = ast.sig;
    let block = ast.block;
    let attrs = ast.attrs;

    let arithmetic: syn::Block = syn::parse(define_arithmetic()).unwrap();
    let intrinsics: syn::Block = syn::parse(define_intrinsics()).unwrap();

    let begin: syn::Stmt = syn::parse((quote! { begin(); }).into()).unwrap();
    let end: syn::Stmt = syn::parse((quote! { end(); }).into()).unwrap();
    let ret: syn::Stmt =
        syn::parse((quote! { return llvm_asm_string.into_inner(); }).into()).unwrap();

    let mut new_stmts = Vec::new();
    for stmt in &intrinsics.stmts {
        new_stmts.push(stmt.clone());
    }
    for stmt in &arithmetic.stmts {
        new_stmts.push(stmt.clone());
    }

    new_stmts.push(begin);

    for stmt in block.stmts {
        new_stmts.push(stmt);
    }

    new_stmts.push(end);
    new_stmts.push(ret);

    let new_block = syn::Block {
        brace_token: block.brace_token,
        stmts: new_stmts,
    };

    let gen = quote! {
        #(#attrs)
        *
        #sig {
            let mut llvm_asm_string = RefCell::new(String::new());

            #new_block
        }
    };
    gen.into()
}
