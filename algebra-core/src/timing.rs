#[macro_export]
macro_rules! timer_println {
    ($now: ident, $string: expr) => {
        #[cfg(feature = "timing")]
        {
            println!("[{:^24}] {} us", $string, $now.1.elapsed().as_micros(),);
        }

        #[cfg(feature = "timing_detailed")]
        {
            macro_rules! function {
                () => {{
                    fn f() {}
                    fn type_name_of<T>(_: T) -> &'static str {
                        core::any::type_name::<T>()
                    }
                    let name = type_name_of(f);
                    &name[..name.len() - 3]
                }};
            }
            let func_string = String::from(function!());
            let mut func_str_vec: Vec<_> = func_string.split("::").collect();
            while *func_str_vec.last().unwrap() == "{{closure}}" {
                func_str_vec.pop();
            }
            println!(
                "{:30} {:26} [{:^24}] {} us",
                format!(
                    "{} {}:{}",
                    String::from(file!()).split("/").last().unwrap(),
                    $now.0,
                    line!(),
                ),
                func_str_vec.last().unwrap(),
                $string,
                $now.1.elapsed().as_micros(),
            );
        }
    };
}

#[macro_export]
macro_rules! timer {
    () => {{
        #[cfg(any(feature = "timing", feature = "timing_detailed"))]
        let now = (line!(), std::time::Instant::now());

        #[cfg(not(any(feature = "timing", feature = "timing_detailed")))]
        let now = ();
        now
    }};
}
