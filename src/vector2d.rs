use std::ops::{Add, Div, Mul, Neg, Sub};

#[derive(Debug, Copy, Clone, PartialEq)]
pub struct Vector2d {
    pub x: f32,
    pub y: f32,
}

impl Vector2d {
    pub fn new(x: f32, y: f32) -> Self {
        Vector2d { x, y }
    }

    pub fn zero() -> Self {
        Vector2d { x: 0.0, y: 0.0 }
    }

    pub fn x_axis() -> Self {
        Vector2d { x: 1.0, y: 0.0 }
    }

    pub fn y_axis() -> Self {
        Vector2d { x: 0.0, y: 1.0 }
    }

    pub fn dot(lhs: Vector2d, rhs: Vector2d) -> f32 {
        (lhs.x * rhs.x) + (lhs.y * rhs.y)
    }

    pub fn distance_squared(lhs: Vector2d, rhs: Vector2d) -> f32 {
        let x_diff = lhs.x - rhs.x;
        let y_diff = lhs.y - rhs.y;

        (x_diff * x_diff) + (y_diff * y_diff)
    }

    pub fn distance(lhs: Vector2d, rhs: Vector2d) -> f32 {
        f32::sqrt(Vector2d::distance_squared(lhs, rhs))
    }

    pub fn angle_in_radians(lhs: Vector2d, rhs: Vector2d) -> f32 {
        let dot_product = Vector2d::dot(lhs, rhs);

        let lhs_magnitude = lhs.magnitude();
        let rhs_magntiude = rhs.magnitude();

        f32::acos(dot_product / (lhs_magnitude * rhs_magntiude))
    }

    pub fn angle_in_degrees(lhs: Vector2d, rhs: Vector2d) -> f32 {
        f32::to_degrees(Vector2d::angle_in_radians(lhs, rhs))
    }

    pub fn unit_vector(v: Vector2d) -> Vector2d {
        let magnitude = v.magnitude();

        let mut unit_vector = v;
        unit_vector.x /= magnitude;
        unit_vector.y /= magnitude;

        unit_vector
    }

    pub fn magnitude(self) -> f32 {
        f32::sqrt(self.magnitude_squared())
    }

    pub fn magnitude_squared(self) -> f32 {
        Vector2d::dot(self, self)
    }

    pub fn normalize(&mut self) {
        let magnitude = self.magnitude();

        self.x /= magnitude;
        self.y /= magnitude;
    }
}

impl Add<Vector2d> for Vector2d {
    type Output = Vector2d;

    fn add(self, rhs: Vector2d) -> Vector2d {
        Vector2d {
            x: self.x + rhs.x,
            y: self.y + rhs.y,
        }
    }
}

impl Sub<Vector2d> for Vector2d {
    type Output = Vector2d;

    fn sub(self, rhs: Vector2d) -> Vector2d {
        Vector2d {
            x: self.x - rhs.x,
            y: self.y - rhs.y,
        }
    }
}

impl Mul<f32> for Vector2d {
    type Output = Vector2d;

    fn mul(self, rhs: f32) -> Vector2d {
        Vector2d {
            x: self.x * rhs,
            y: self.y * rhs,
        }
    }
}

impl Div<f32> for Vector2d {
    type Output = Vector2d;

    fn div(self, rhs: f32) -> Vector2d {
        if rhs == 0.0 {
            panic!("Zero is an invalid denominator!")
        }

        Vector2d {
            x: self.x / rhs,
            y: self.y / rhs,
        }
    }
}

impl Neg for Vector2d {
    type Output = Vector2d;

    fn neg(self) -> Vector2d {
        Vector2d {
            x: -self.x,
            y: -self.y,
        }
    }
}
