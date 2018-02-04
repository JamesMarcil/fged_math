use std::ops::{Add, Div, Mul, Neg, Sub};

#[derive(Debug, Copy, Clone, PartialEq)]
pub struct Vector4d {
    pub x: f32,
    pub y: f32,
    pub z: f32,
    pub w: f32,
}

impl Vector4d {
    pub fn new(x: f32, y: f32, z: f32, w:f32) -> Self {
        Vector4d { x, y, z, w}
    }

    pub fn zero() -> Self {
        Vector4d {
            x: 0.0,
            y: 0.0,
            z: 0.0,
            w: 0.0,
        }
    }

    pub fn x_axis() -> Self {
        Vector4d {
            x: 1.0,
            y: 0.0,
            z: 0.0,
            w: 0.0,
        }
    }

    pub fn y_axis() -> Self {
        Vector4d {
            x: 0.0,
            y: 1.0,
            z: 0.0,
            w: 0.0,
        }
    }

    pub fn z_axis() -> Self {
        Vector4d {
            x: 0.0,
            y: 0.0,
            z: 1.0,
            w: 0.0,
        }
    }

    pub fn w_axis() -> Self {
        Vector4d {
            x: 0.0,
            y: 0.0,
            z: 0.0,
            w: 1.0,
        }
    }

    pub fn dot(lhs: Vector4d, rhs: Vector4d) -> f32 {
        (lhs.x * rhs.x) + (lhs.y * rhs.y) + (lhs.z * rhs.z) + (lhs.w * rhs.w)
    }

    pub fn distance_squared(lhs: Vector4d, rhs: Vector4d) -> f32 {
        let x_diff = lhs.x - rhs.x;
        let y_diff = lhs.y - rhs.y;
        let z_diff = lhs.z - rhs.z;
        let w_diff = lhs.w - rhs.w;

        (x_diff * x_diff) + (y_diff * y_diff) + (z_diff * z_diff) + (w_diff * w_diff)
    }

    pub fn distance(lhs: Vector4d, rhs: Vector4d) -> f32 {
        f32::sqrt(Vector4d::distance_squared(lhs, rhs))
    }

    pub fn angle_in_radians(lhs: Vector4d, rhs: Vector4d) -> f32 {
        let dot_product = Vector4d::dot(lhs, rhs);

        let lhs_magnitude = lhs.magnitude();
        let rhs_magntiude = rhs.magnitude();

        f32::acos(dot_product / (lhs_magnitude * rhs_magntiude))
    }

    pub fn angle_in_degrees(lhs: Vector4d, rhs: Vector4d) -> f32 {
        f32::to_degrees(Vector4d::angle_in_radians(lhs, rhs))
    }

    pub fn unit_vector(v: Vector4d) -> Vector4d {
        let magnitude = v.magnitude();

        let mut unit_vector = v;
        unit_vector.x /= magnitude;
        unit_vector.y /= magnitude;
        unit_vector.z /= magnitude;
        unit_vector.w /= magnitude;

        unit_vector
    }

    pub fn magnitude(self) -> f32 {
        f32::sqrt(self.magnitude_squared())
    }

    pub fn magnitude_squared(self) -> f32 {
        Vector4d::dot(self, self)
    }

    pub fn normalize(&mut self) {
        let magnitude = self.magnitude();

        self.x /= magnitude;
        self.y /= magnitude;
        self.z /= magnitude;
        self.w /= magnitude;
    }
}

impl Add<Vector4d> for Vector4d {
    type Output = Vector4d;

    fn add(self, rhs: Vector4d) -> Vector4d {
        Vector4d {
            x: self.x + rhs.x,
            y: self.y + rhs.y,
            z: self.z + rhs.z,
            w: self.w + rhs.w,
        }
    }
}

impl Sub<Vector4d> for Vector4d {
    type Output = Vector4d;

    fn sub(self, rhs: Vector4d) -> Vector4d {
        Vector4d {
            x: self.x - rhs.x,
            y: self.y - rhs.y,
            z: self.z - rhs.z,
            w: self.w - rhs.w,
        }
    }
}

impl Mul<f32> for Vector4d {
    type Output = Vector4d;

    fn mul(self, rhs: f32) -> Vector4d {
        Vector4d {
            x: self.x * rhs,
            y: self.y * rhs,
            z: self.z * rhs,
            w: self.w * rhs,
        }
    }
}

impl Div<f32> for Vector4d {
    type Output = Vector4d;

    fn div(self, rhs: f32) -> Vector4d {
        Vector4d {
            x: self.x / rhs,
            y: self.y / rhs,
            z: self.z / rhs,
            w: self.w / rhs,
        }
    }
}

impl Neg for Vector4d {
    type Output = Vector4d;

    fn neg(self) -> Vector4d {
        Vector4d {
            x: -self.x,
            y: -self.y,
            z: -self.z,
            w: -self.w,
        }
    }
}
