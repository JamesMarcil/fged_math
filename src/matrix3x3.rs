use std::ops::{Add, Index, IndexMut, Mul, Div, Sub};

#[derive(Debug, Clone, Copy)]
pub struct Matrix3x3 {
    elements: [f32; 9],
}

impl Matrix3x3 {
    pub fn zero() -> Matrix3x3 {
        Matrix3x3 { elements: [0.0; 9] }
    }

    pub fn new(
        m00: f32,
        m01: f32,
        m02: f32,
        m10: f32,
        m11: f32,
        m12: f32,
        m20: f32,
        m21: f32,
        m22: f32,
    ) -> Matrix3x3 {
        let mut matrix = Matrix3x3::zero();

        matrix[0] = m00;
        matrix[1] = m10;
        matrix[2] = m20;

        matrix[3] = m01;
        matrix[4] = m11;
        matrix[5] = m21;

        matrix[6] = m02;
        matrix[7] = m12;
        matrix[8] = m22;

        matrix
    }

    pub fn identity() -> Matrix3x3 {
        let mut identity = Matrix3x3::zero();

        identity[0] = 1.0;
        identity[4] = 1.0;
        identity[8] = 1.0;

        identity
    }
}

impl Add<Matrix3x3> for Matrix3x3 {
    type Output = Matrix3x3;

    fn add(self, rhs: Matrix3x3) -> Matrix3x3 {
        let mut matrix = Matrix3x3::zero();

        matrix[0] = self[0] + rhs[0];
        matrix[1] = self[1] + rhs[1];
        matrix[2] = self[2] + rhs[2];
        matrix[3] = self[3] + rhs[3];
        matrix[4] = self[4] + rhs[4];
        matrix[5] = self[5] + rhs[5];
        matrix[6] = self[6] + rhs[6];
        matrix[7] = self[7] + rhs[7];
        matrix[8] = self[8] + rhs[8];

        matrix
    }
}

impl Sub<Matrix3x3> for Matrix3x3 {
    type Output = Matrix3x3;

    fn sub(self, rhs: Matrix3x3) -> Matrix3x3 {
        let mut matrix = Matrix3x3::zero();

        matrix[0] = self[0] - rhs[0];
        matrix[1] = self[1] - rhs[1];
        matrix[2] = self[2] - rhs[2];
        matrix[3] = self[3] - rhs[3];
        matrix[4] = self[4] - rhs[4];
        matrix[5] = self[5] - rhs[5];
        matrix[6] = self[6] - rhs[6];
        matrix[7] = self[7] - rhs[7];
        matrix[8] = self[8] - rhs[8];

        matrix
    }
}

impl Mul<f32> for Matrix3x3 {
    type Output = Matrix3x3;

    fn mul(self, rhs: f32) -> Matrix3x3 {
        let mut matrix = Matrix3x3::zero();

        matrix[0] = self[0] * rhs;
        matrix[1] = self[1] * rhs;
        matrix[2] = self[2] * rhs;
        matrix[3] = self[3] * rhs;
        matrix[4] = self[4] * rhs;
        matrix[5] = self[5] * rhs;
        matrix[6] = self[6] * rhs;
        matrix[7] = self[7] * rhs;
        matrix[8] = self[8] * rhs;

        matrix
    }
}

impl Div<f32> for Matrix3x3 {
    type Output = Matrix3x3;

    fn div(self, rhs: f32) -> Matrix3x3 {
        let mut matrix = Matrix3x3::zero();

        matrix[0] = self[0] / rhs;
        matrix[1] = self[1] / rhs;
        matrix[2] = self[2] / rhs;
        matrix[3] = self[3] / rhs;
        matrix[4] = self[4] / rhs;
        matrix[5] = self[5] / rhs;
        matrix[6] = self[6] / rhs;
        matrix[7] = self[7] / rhs;
        matrix[8] = self[8] / rhs;

        matrix
    }
}

impl Mul<Matrix3x3> for Matrix3x3 {
    type Output = Matrix3x3;

    fn mul(self, rhs: Matrix3x3) -> Matrix3x3 {
        let mut matrix = Matrix3x3::zero();

        matrix[0] = self[0] * rhs[0] + self[3] * rhs[1] + self[6] * rhs[2]; // m00 = a00 * b00 + a01 * b10 + a02 * b20
        matrix[1] = self[1] * rhs[0] + self[4] * rhs[1] + self[7] * rhs[2]; // m10 = a10 * b00 + a11 * b10 + a12 * b20
        matrix[2] = self[2] * rhs[0] + self[5] * rhs[1] + self[8] * rhs[2]; // m20 = a20 * b00 + a21 * b10 + a22 * b20

        matrix[3] = self[0] * rhs[3] + self[3] * rhs[4] + self[6] * rhs[5]; // m01 = a00 * b01 + a01 * b11 + a02 * b21
        matrix[4] = self[1] * rhs[3] + self[4] * rhs[4] + self[7] * rhs[5]; // m11 = a10 * b01 + a11 * b11 + a12 * b21
        matrix[5] = self[2] * rhs[3] + self[5] * rhs[4] + self[8] * rhs[5]; // m12 = a20 * b01 + a21 * b11 + a22 * b21

        matrix[6] = self[0] * rhs[6] + self[3] * rhs[7] + self[6] * rhs[8]; // m02 = a00 * b02 + a01 * b12 + a02 * b22
        matrix[7] = self[1] * rhs[6] + self[4] * rhs[7] + self[7] * rhs[8]; // m21 = a10 * b02 + a11 * b12 + a12 * b22
        matrix[8] = self[2] * rhs[6] + self[5] * rhs[7] + self[8] * rhs[8]; // m22 = a20 * b02 + a21 * b12 + a22 * b22

        matrix
    }
}

impl Index<usize> for Matrix3x3 {
    type Output = f32;

    fn index(&self, index: usize) -> &f32 {
        &self.elements[index]
    }
}

impl IndexMut<usize> for Matrix3x3 {
    fn index_mut(&mut self, index: usize) -> &mut f32 {
        &mut self.elements[index]
    }
}
