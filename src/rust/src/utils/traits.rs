use bincode::{Decode, Encode};
use extendr_api::*;
use faer_entity::SimpleEntity;
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

impl std::ops::Add for F16 {
    type Output = F16;

    fn add(self, other: F16) -> F16 {
        let a = f16::from(self);
        let b = f16::from(other);
        F16::from(a + b)
    }
}

impl<'a> std::ops::Add<&'a F16> for F16 {
    type Output = F16;

    fn add(self, other: &'a F16) -> F16 {
        let a = f16::from(self);
        let b = f16::from(*other);
        F16::from(a + b)
    }
}

impl std::ops::Mul for F16 {
    type Output = F16;

    fn mul(self, other: F16) -> F16 {
        let a = f16::from(self);
        let b = f16::from(other);
        F16::from(a * b)
    }
}

impl<'a> std::ops::Mul<&'a F16> for F16 {
    type Output = F16;

    fn mul(self, other: &'a F16) -> F16 {
        let a = f16::from(self);
        let b = f16::from(*other);
        F16::from(a * b)
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

    /// Create an F16
    pub fn to_f32(self) -> f32 {
        f16::from_bits(self.0).to_f32()
    }
}

///////////////////////
// f32 and u16 stuff //
///////////////////////

pub trait ToF32AndU16: Copy {
    fn to_f32(self) -> f32;
    fn to_u16(self) -> u16;
}

impl ToF32AndU16 for i32 {
    fn to_f32(self) -> f32 {
        self as f32
    }
    fn to_u16(self) -> u16 {
        self as u16
    }
}

impl ToF32AndU16 for u32 {
    fn to_f32(self) -> f32 {
        self as f32
    }
    fn to_u16(self) -> u16 {
        self as u16
    }
}

impl ToF32AndU16 for usize {
    fn to_f32(self) -> f32 {
        self as f32
    }
    fn to_u16(self) -> u16 {
        self as u16
    }
}

impl ToF32AndU16 for u16 {
    fn to_f32(self) -> f32 {
        self as f32
    }
    fn to_u16(self) -> u16 {
        self
    }
}

////////////////////
// R vector stuff //
////////////////////

pub trait VecConvert<U> {
    fn r_int_convert(self) -> Vec<U>;
}

impl VecConvert<i32> for Vec<usize> {
    fn r_int_convert(self) -> Vec<i32> {
        self.into_iter().map(|x| x as i32).collect()
    }
}

impl VecConvert<i32> for &[&usize] {
    fn r_int_convert(self) -> Vec<i32> {
        self.iter().map(|&&x| x as i32).collect()
    }
}

impl VecConvert<i32> for Vec<u32> {
    fn r_int_convert(self) -> Vec<i32> {
        self.into_iter().map(|x| x as i32).collect()
    }
}

impl VecConvert<usize> for Vec<i32> {
    fn r_int_convert(self) -> Vec<usize> {
        self.into_iter().map(|x| x as usize).collect()
    }
}

pub trait VecFloatConvert<U> {
    fn r_float_convert(self) -> Vec<U>;
}

impl VecFloatConvert<f32> for Vec<f64> {
    fn r_float_convert(self) -> Vec<f32> {
        self.into_iter().map(|x| x as f32).collect()
    }
}

impl VecFloatConvert<f64> for Vec<f32> {
    fn r_float_convert(self) -> Vec<f64> {
        self.into_iter().map(|x| x as f64).collect()
    }
}

impl VecFloatConvert<f64> for &[f32] {
    fn r_float_convert(self) -> Vec<f64> {
        self.iter().map(|x| *x as f64).collect()
    }
}

////////////////
// R and faer //
////////////////

pub trait FaerRType: SimpleEntity + Copy + Clone + 'static {
    type RType: Copy + Clone;
    fn to_r_matrix(x: faer::MatRef<Self>) -> extendr_api::RArray<Self::RType, [usize; 2]>;
}

impl FaerRType for f64 {
    type RType = f64;
    fn to_r_matrix(x: faer::MatRef<Self>) -> extendr_api::RArray<Self, [usize; 2]> {
        let nrow = x.nrows();
        let ncol = x.ncols();
        RArray::new_matrix(nrow, ncol, |row, column| x[(row, column)])
    }
}

impl FaerRType for i32 {
    type RType = i32;
    fn to_r_matrix(x: faer::MatRef<Self>) -> extendr_api::RArray<Self, [usize; 2]> {
        let nrow = x.nrows();
        let ncol = x.ncols();
        RArray::new_matrix(nrow, ncol, |row, column| x[(row, column)])
    }
}

impl FaerRType for f32 {
    type RType = f64;
    fn to_r_matrix(x: faer::MatRef<Self>) -> extendr_api::RArray<f64, [usize; 2]> {
        let nrow = x.nrows();
        let ncol = x.ncols();
        RArray::new_matrix(nrow, ncol, |row, column| x[(row, column)] as f64)
    }
}
