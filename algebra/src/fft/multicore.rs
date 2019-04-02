//! This is an interface for dealing with the kinds of
//! parallel computations involved in `snark`. It's
//! currently just a thin wrapper around `rayon`.
use rayon::{self, Scope};

#[derive(Clone)]
pub struct Worker {
    cpus: usize,
}

impl Worker {
    pub fn new() -> Worker {
        let cpus = rayon::current_num_threads();
        Self { cpus }
    }

    pub fn log_num_cpus(&self) -> u32 {
        log2_floor(self.cpus)
    }

    pub fn scope<'a, F, R>(&self, elements: usize, f: F) -> R
    where
        F: 'a + Send + FnOnce(&Scope<'a>, usize) -> R,
        R: Send,
    {
        let chunk_size = if elements < self.cpus {
            1
        } else {
            elements / self.cpus
        };

        rayon::scope(move |scope| f(scope, chunk_size))
    }
}

fn log2_floor(num: usize) -> u32 {
    assert!(num > 0);

    let mut pow = 0;

    while (1 << (pow + 1)) <= num {
        pow += 1;
    }

    pow
}

#[test]
fn test_log2_floor() {
    assert_eq!(log2_floor(1), 0);
    assert_eq!(log2_floor(2), 1);
    assert_eq!(log2_floor(3), 1);
    assert_eq!(log2_floor(4), 2);
    assert_eq!(log2_floor(5), 2);
    assert_eq!(log2_floor(6), 2);
    assert_eq!(log2_floor(7), 2);
    assert_eq!(log2_floor(8), 3);
}
