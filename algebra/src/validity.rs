pub trait SemanticallyValid {
    /// Does all the necessary checks on `Self` to estabilish its semantic validity.
    /// NOTE: The meaning of "semantic validity" for `Self` is actually defined by
    /// the implementation of this function.
    fn is_valid(&self) -> bool;
}

impl<T: SemanticallyValid> SemanticallyValid for Vec<T> {
    fn is_valid(&self) -> bool {
        for item in self.iter() {
            if !item.is_valid() {
                return false;
            }
        }
        true
    }
}
