//! no-std io replacement
use crate::Vec;
use core::{cmp, fmt, mem};

#[derive(Debug)]
pub struct Error;

impl fmt::Display for Error {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        writeln!(f, "Error")
    }
}

pub type Result<T> = core::result::Result<T, Error>;

pub trait Read {
    fn read_exact(&mut self, data: &mut [u8]) -> Result<()>;
}

pub trait Write {
    fn write_all(&mut self, data: &[u8]) -> Result<()>;
}

impl<R: Read + ?Sized> Read for &mut R {
    #[inline]
    fn read_exact(&mut self, data: &mut [u8]) -> Result<()> {
        (**self).read_exact(data)
    }
}

impl Read for &[u8] {
    fn read_exact(&mut self, buf: &mut [u8]) -> Result<()> {
        if buf.len() > self.len() {
            return Err(Error);
        }
        let (a, b) = self.split_at(buf.len());

        // First check if the amount of bytes we want to read is small:
        // `copy_from_slice` will generally expand to a call to `memcpy`, and
        // for a single byte the overhead is significant.
        if buf.len() == 1 {
            buf[0] = a[0];
        } else {
            buf.copy_from_slice(a);
        }

        *self = b;
        Ok(())
    }
}

impl<W: Write + ?Sized> Write for &mut W {
    #[inline]
    fn write_all(&mut self, data: &[u8]) -> Result<()> {
        (**self).write_all(data)
    }
}

impl Write for &mut [u8] {
    #[inline]
    fn write_all(&mut self, data: &[u8]) -> Result<()> {
        let amt = cmp::min(data.len(), self.len());
        let (a, b) = mem::replace(self, &mut []).split_at_mut(amt);
        a.copy_from_slice(&data[..amt]);
        *self = b;

        if amt == data.len() {
            Ok(())
        } else {
            Err(Error)
        }
    }
}

impl Write for Vec<u8> {
    #[inline]
    fn write_all(&mut self, buf: &[u8]) -> Result<()> {
        self.extend_from_slice(buf);
        Ok(())
    }
}

/// This data structure is used as a workaround for current design of `ToBytes`
/// which does not allow multiple writes to `&mut [u8]`.
pub struct Cursor<T> {
    inner: T,
    pos: usize,
}

impl<T> Cursor<T> {
    pub fn new(inner: T) -> Self {
        Cursor { inner, pos: 0 }
    }
}

impl Write for Cursor<&mut [u8]> {
    fn write_all(&mut self, buf: &[u8]) -> Result<()> {
        let to_copy = cmp::min(self.inner.len() - self.pos, buf.len());
        self.inner[self.pos..self.pos + to_copy].copy_from_slice(buf);
        self.pos += to_copy;
        Ok(())
    }
}

impl Read for Cursor<&mut [u8]> {
    fn read_exact(&mut self, buf: &mut [u8]) -> Result<()> {
        let to_copy = cmp::min(self.inner.len() - self.pos, buf.len());
        buf.copy_from_slice(&self.inner[self.pos..self.pos + to_copy]);
        self.pos += to_copy;
        Ok(())
    }
}

impl Read for Cursor<&[u8]> {
    fn read_exact(&mut self, buf: &mut [u8]) -> Result<()> {
        let to_copy = cmp::min(self.inner.len() - self.pos, buf.len());
        buf.copy_from_slice(&self.inner[self.pos..self.pos + to_copy]);
        self.pos += to_copy;
        Ok(())
    }
}
