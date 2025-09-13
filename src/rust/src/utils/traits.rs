use bincode::{Decode, Encode};
use half::f16;
use serde::{Deserialize, Serialize};
use std::cmp::Ordering;
use std::iter::Sum;

/////////
// F16 //
/////////

#[derive(Encode, Decode, Serialize, Deserialize, Debug, Clone, Copy, Default)]
pub struct F16(u16);

impl From<f16> for F16 {
    fn from(f: f16) -> Self {
        F16(f.to_bits())
    }
}

impl From<F16> for f16 {
    fn from(f: F16) -> Self {
        f16::from_bits(f.0)
    }
}

impl Sum for F16 {
    fn sum<I: Iterator<Item = F16>>(iter: I) -> Self {
        let sum: f16 = iter.map(f16::from).sum();
        F16::from(sum)
    }
}

impl<'a> Sum<&'a F16> for F16 {
    fn sum<I: Iterator<Item = &'a F16>>(iter: I) -> Self {
        let sum: f16 = iter.map(|&f| f16::from(f)).sum();
        F16::from(sum)
    }
}

impl PartialEq for F16 {
    fn eq(&self, other: &Self) -> bool {
        let a = f16::from(*self);
        let b = f16::from(*other);

        if a.is_nan() && b.is_nan() {
            false
        } else {
            a == b
        }
    }
}

impl Eq for F16 {}

impl PartialOrd for F16 {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for F16 {
    fn cmp(&self, other: &Self) -> Ordering {
        let a = half::f16::from(*self);
        let b = half::f16::from(*other);

        match (a.is_nan(), b.is_nan()) {
            (true, true) => Ordering::Equal,
            (true, false) => Ordering::Greater,
            (false, true) => Ordering::Less,
            (false, false) => a.partial_cmp(&b).unwrap(),
        }
    }
}

impl F16 {
    /// Get the raw bits representation
    pub fn to_bits(self) -> u16 {
        self.0
    }

    /// Create from raw bits
    pub fn from_bits(bits: u16) -> Self {
        F16(bits)
    }

    /// From f32
    pub fn from_f32(f: f32) -> Self {
        F16::from(f16::from_f32(f))
    }

    /// Convert to little-endian bytes
    pub fn to_le_bytes(self) -> [u8; 2] {
        self.0.to_le_bytes()
    }

    /// Create from little-endian bytes
    pub fn from_le_bytes(bytes: [u8; 2]) -> Self {
        F16(u16::from_le_bytes(bytes))
    }
}
