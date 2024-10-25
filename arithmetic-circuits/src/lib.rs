pub mod arithmetic_circuit;
pub mod expression;
pub mod reader;

#[macro_export]
macro_rules! TEST_DATA_PATH {
    () => {
        concat!(env!("CARGO_MANIFEST_DIR"), "/../circom/{}",)
    };
}
