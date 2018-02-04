use std::ops::{Add, Div, Mul, Neg, Sub};

#[derive(Debug, Copy, Clone, PartialEq)]
pub struct Vector3d {
    pub x: f32,
    pub y: f32,
    pub z: f32,
}

impl Vector3d {
    pub fn new(x: f32, y: f32, z: f32) -> Self {
        Vector3d { x, y, z }
    }

    pub fn zero() -> Self {
        Vector3d {
            x: 0.0,
            y: 0.0,
            z: 0.0,
        }
    }

    pub fn x_axis() -> Self {
        Vector3d {
            x: 1.0,
            y: 0.0,
            z: 0.0,
        }
    }

    pub fn y_axis() -> Self {
        Vector3d {
            x: 0.0,
            y: 1.0,
            z: 0.0,
        }
    }

    pub fn z_axis() -> Self {
        Vector3d {
            x: 0.0,
            y: 0.0,
            z: 1.0,
        }
    }

    pub fn dot(lhs: Vector3d, rhs: Vector3d) -> f32 {
        (lhs.x * rhs.x) + (lhs.y * rhs.y) + (lhs.z * rhs.z)
    }

    pub fn distance_squared(lhs: Vector3d, rhs: Vector3d) -> f32 {
        let x_diff = lhs.x - rhs.x;
        let y_diff = lhs.y - rhs.y;
        let z_diff = lhs.z - rhs.z;

        (x_diff * x_diff) + (y_diff * y_diff) + (z_diff * z_diff)
    }

    pub fn distance(lhs: Vector3d, rhs: Vector3d) -> f32 {
        f32::sqrt(Vector3d::distance_squared(lhs, rhs))
    }

    pub fn angle_in_radians(lhs: Vector3d, rhs: Vector3d) -> f32 {
        let dot_product = Vector3d::dot(lhs, rhs);

        let lhs_magnitude = lhs.magnitude();
        let rhs_magntiude = rhs.magnitude();

        f32::acos(dot_product / (lhs_magnitude * rhs_magntiude))
    }

    pub fn angle_in_degrees(lhs: Vector3d, rhs: Vector3d) -> f32 {
        f32::to_degrees(Vector3d::angle_in_radians(lhs, rhs))
    }

    pub fn unit_vector(v: Vector3d) -> Vector3d {
        let magnitude = v.magnitude();

        let mut unit_vector = v;
        unit_vector.x /= magnitude;
        unit_vector.y /= magnitude;
        unit_vector.z /= magnitude;

        unit_vector
    }

    pub fn magnitude(self) -> f32 {
        f32::sqrt(self.magnitude_squared())
    }

    pub fn magnitude_squared(self) -> f32 {
        Vector3d::dot(self, self)
    }

    pub fn normalize(&mut self) {
        let magnitude = self.magnitude();

        self.x /= magnitude;
        self.y /= magnitude;
        self.z /= magnitude;
    }
}

impl Add<Vector3d> for Vector3d {
    type Output = Vector3d;

    fn add(self, rhs: Vector3d) -> Vector3d {
        Vector3d {
            x: self.x + rhs.x,
            y: self.y + rhs.y,
            z: self.z + rhs.z,
        }
    }
}

impl Sub<Vector3d> for Vector3d {
    type Output = Vector3d;

    fn sub(self, rhs: Vector3d) -> Vector3d {
        Vector3d {
            x: self.x - rhs.x,
            y: self.y - rhs.y,
            z: self.z - rhs.z,
        }
    }
}

impl Mul<f32> for Vector3d {
    type Output = Vector3d;

    fn mul(self, rhs: f32) -> Vector3d {
        Vector3d {
            x: self.x * rhs,
            y: self.y * rhs,
            z: self.z * rhs,
        }
    }
}

impl Div<f32> for Vector3d {
    type Output = Vector3d;

    fn div(self, rhs: f32) -> Vector3d {
        Vector3d {
            x: self.x / rhs,
            y: self.y / rhs,
            z: self.z / rhs,
        }
    }
}

impl Neg for Vector3d {
    type Output = Vector3d;

    fn neg(self) -> Vector3d {
        Vector3d {
            x: -self.x,
            y: -self.y,
            z: -self.z,
        }
    }
}
