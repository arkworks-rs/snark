// adapted from `tracing_error::{SpanTrace, ErrorLayer}`.

use core::any::{type_name, TypeId};
use core::fmt;
use core::marker::PhantomData;
use tracing::{span, Dispatch, Metadata, Subscriber};
use tracing_subscriber::{
    layer::{self, Layer},
    registry::LookupSpan,
};

/// A subscriber [`Layer`] that enables capturing a trace of R1CS constraint generation.
///
/// [`Layer`]: https://docs.rs/tracing-subscriber/0.2.10/tracing_subscriber/layer/trait.Layer.html
/// [field formatter]: https://docs.rs/tracing-subscriber/0.2.10/tracing_subscriber/fmt/trait.FormatFields.html
/// [default format]: https://docs.rs/tracing-subscriber/0.2.10/tracing_subscriber/fmt/format/struct.DefaultFields.html
pub struct ConstraintLayer<S> {
    /// Mode of filtering.
    pub mode: TracingMode,

    get_context: WithContext,
    _subscriber: PhantomData<fn(S)>,
}

/// Instructs `ConstraintLayer` to conditionally filter out spans.
#[derive(PartialEq, Eq, Ord, PartialOrd, Hash, Debug)]
pub enum TracingMode {
    /// Instructs `ConstraintLayer` to filter out any spans that *do not* have `target == "r1cs"`.
    OnlyConstraints,
    /// Instructs `ConstraintLayer` to filter out any spans that *do* have `target == "r1cs"`.
    NoConstraints,
    /// Instructs `ConstraintLayer` to not filter out any spans.
    All,
}

// this function "remembers" the types of the subscriber and the formatter,
// so that we can downcast to something aware of them without knowing those
// types at the callsite.
pub(crate) struct WithContext(
    fn(&Dispatch, &span::Id, f: &mut dyn FnMut(&'static Metadata<'static>, &str) -> bool),
);

impl<S> Layer<S> for ConstraintLayer<S>
where
    S: Subscriber + for<'span> LookupSpan<'span>,
{
    fn enabled(&self, metadata: &Metadata, _ctx: layer::Context<S>) -> bool {
        match self.mode {
            TracingMode::OnlyConstraints => metadata.target() == "r1cs",
            TracingMode::NoConstraints => metadata.target() != "r1cs",
            TracingMode::All => true,
        }
    }

    /// Notifies this layer that a new span was constructed with the given
    /// `Attributes` and `Id`.
    fn new_span(&self, _attrs: &span::Attributes<'_>, _id: &span::Id, _ctx: layer::Context<'_, S>) {
    }

    #[allow(unsafe_code, trivial_casts)]
    unsafe fn downcast_raw(&self, id: TypeId) -> Option<*const ()> {
        match id {
            id if id == TypeId::of::<Self>() => Some(self as *const _ as *const ()),
            id if id == TypeId::of::<WithContext>() => {
                Some(&self.get_context as *const _ as *const ())
            }
            _ => None,
        }
    }
}

impl<S> ConstraintLayer<S>
where
    S: Subscriber + for<'span> LookupSpan<'span>,
{
    /// Returns a new `ConstraintLayer`.
    ///
    /// If `mode == TracingMode::OnlyConstraints`, the resulting layer will
    /// filter out any spans whose `target != "r1cs"`.
    ///
    /// If `mode == TracingMode::NoConstraints`, the resulting layer will
    /// filter out any spans whose `target == "r1cs"`.
    ///
    /// Finally, if `mode == TracingMode::All`, the resulting layer will
    /// not filter out any spans.
    pub fn new(mode: TracingMode) -> Self {
        Self {
            mode,
            get_context: WithContext(Self::get_context),
            _subscriber: PhantomData,
        }
    }

    fn get_context(
        dispatch: &Dispatch,
        id: &span::Id,
        f: &mut dyn FnMut(&'static Metadata<'static>, &str) -> bool,
    ) {
        let subscriber = dispatch
            .downcast_ref::<S>()
            .expect("subscriber should downcast to expected type; this is a bug!");
        let span = subscriber
            .span(id)
            .expect("registry should have a span for the current ID");
        let parents = span.parents();
        for span in std::iter::once(span).chain(parents) {
            let cont = f(span.metadata(), "");
            if !cont {
                break;
            }
        }
    }
}

impl WithContext {
    pub(crate) fn with_context<'a>(
        &self,
        dispatch: &'a Dispatch,
        id: &span::Id,
        mut f: impl FnMut(&'static Metadata<'static>, &str) -> bool,
    ) {
        (self.0)(dispatch, id, &mut f)
    }
}

impl<S> Default for ConstraintLayer<S>
where
    S: Subscriber + for<'span> LookupSpan<'span>,
{
    fn default() -> Self {
        Self::new(TracingMode::All)
    }
}

impl<S> fmt::Debug for ConstraintLayer<S> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("ConstraintLayer")
            .field("subscriber", &format_args!("{}", type_name::<S>()))
            .finish()
    }
}

macro_rules! try_bool {
    ($e:expr, $dest:ident) => {{
        let ret = $e.unwrap_or_else(|e| $dest = Err(e));

        if $dest.is_err() {
            return false;
        }

        ret
    }};
}

/// A captured trace of [`tracing`] spans that have `target = "r1cs"`.
///
/// This type can be thought of as a relative of
/// [`std::backtrace::Backtrace`][`Backtrace`].
/// However, rather than capturing the current call stack when it is
/// constructed, a `ConstraintTrace` instead captures the current [span] and its
/// [parents]. It allows inspection of the constraints that are left unsatisfied
/// by a particular witness assignment to an R1CS instance.
///
/// # Formatting
///
/// The `ConstraintTrace` type implements `fmt::Display`, formatting the span trace
/// similarly to how Rust formats panics. For example:
///
/// ```text
///    0: r1cs-std::bits::something
///           at r1cs-std/src/bits/test.rs:42
///    1: r1cs-std::bits::another_thing
///           at r1cs-std/src/bits/test.rs:15
/// ```
///
/// [`tracing`]: https://docs.rs/tracing
/// [`Backtrace`]: https://doc.rust-lang.org/std/backtrace/struct.Backtrace.html
/// [span]: https://docs.rs/tracing/latest/tracing/span/index.html
/// [parents]: https://docs.rs/tracing/latest/tracing/span/index.html#span-relationships
#[derive(Clone, Debug)]
pub struct ConstraintTrace {
    span: span::Span,
}

// === impl ConstraintTrace ===

impl ConstraintTrace {
    /// Capture the current span trace.
    ///
    /// # Examples
    /// ```rust
    /// use r1cs_core::ConstraintTrace;
    ///
    /// pub struct MyError {
    ///     trace: Option<ConstraintTrace>,
    ///     // ...
    /// }
    ///
    /// # fn some_error_condition() -> bool { true }
    ///
    /// pub fn my_function(arg: &str) -> Result<(), MyError> {
    ///     let _span = tracing::info_span!(target: "r1cs", "In my_function");
    ///     let _guard = _span.enter();
    ///     if some_error_condition() {
    ///         return Err(MyError {
    ///             trace: ConstraintTrace::capture(),
    ///             // ...
    ///         });
    ///     }
    ///
    ///     // ...
    /// #   Ok(())
    /// }
    /// ```
    pub fn capture() -> Option<Self> {
        let span = span::Span::current();

        if span.is_none() {
            None
        } else {
            let trace = Self { span };
            Some(trace)
        }
    }

    /// Apply a function to all captured spans in the trace until it returns
    /// `false`.
    ///
    /// This will call the provided function with a reference to the
    /// [`Metadata`] and a formatted representation of the [fields] of each span
    /// captured in the trace, starting with the span that was current when the
    /// trace was captured. The function may return `true` or `false` to
    /// indicate whether to continue iterating over spans; if it returns
    /// `false`, no additional spans will be visited.
    ///
    /// [fields]: https://docs.rs/tracing/latest/tracing/field/index.html
    /// [`Metadata`]: https://docs.rs/tracing/latest/tracing/struct.Metadata.html
    fn with_spans(&self, f: impl FnMut(&'static Metadata<'static>, &str) -> bool) {
        self.span.with_subscriber(|(id, s)| {
            if let Some(getcx) = s.downcast_ref::<WithContext>() {
                getcx.with_context(s, id, f);
            }
        });
    }

    /// Compute a `Vec` of `TraceStep`s, one for each `Span` on the path from the root
    /// `Span`.
    ///
    /// The output starts from the root of the span tree.
    pub fn path(&self) -> Vec<TraceStep> {
        let mut path = Vec::new();
        self.with_spans(|metadata, _| {
            if metadata.target() == "r1cs" {
                let n = metadata.name();
                let step = metadata
                    .module_path()
                    .map(|m| (n, m))
                    .and_then(|(n, m)| metadata.file().map(|f| (n, m, f)))
                    .and_then(|(n, m, f)| metadata.line().map(|l| (n, m, f, l)));
                if let Some((name, module_path, file, line)) = step {
                    let step = TraceStep {
                        name,
                        module_path,
                        file,
                        line,
                    };
                    path.push(step);
                } else {
                    return false;
                }
            }
            true
        });
        path.reverse(); // root first
        path
    }
}

impl fmt::Display for ConstraintTrace {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let mut err = Ok(());
        let mut span = 0;

        self.with_spans(|metadata, _| {
            if metadata.target() != "r1cs" {
                return true;
            }
            if span > 0 {
                try_bool!(write!(f, "\n",), err);
            }

            try_bool!(
                write!(
                    f,
                    "{:>4}: {}::{}",
                    span,
                    metadata.module_path().unwrap(),
                    metadata.name()
                ),
                err
            );

            if let Some((file, line)) = metadata
                .file()
                .and_then(|file| metadata.line().map(|line| (file, line)))
            {
                try_bool!(write!(f, "\n             at {}:{}", file, line), err);
            }

            span += 1;
            true
        });

        err
    }
}
/// A step in the trace of a constraint generation step.
#[derive(Debug, Clone, Copy)]
pub struct TraceStep {
    /// Name of the constraint generating span.
    pub name: &'static str,
    /// Name of the module containing the constraint generating span.
    pub module_path: &'static str,
    /// Name of the file containing the constraint generating span.
    pub file: &'static str,
    /// Line number of the constraint generating span.
    pub line: u32,
}
