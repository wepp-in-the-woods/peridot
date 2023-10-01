#[allow(dead_code)]
pub fn argmax(data: &Vec<f64>) -> Option<usize> {
    data.iter().enumerate().fold(None, |acc, (index, &value)| {
        match acc {
            None => Some((index, value)),
            Some((_, max_val)) if value > max_val => Some((index, value)),
            _ => acc
        }
    }).map(|(index, _)| index)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_argmax_basic() {
        let data = vec![1.0, 2.0, 5.0, 3.0, 4.0];
        assert_eq!(argmax(&data), Some(2));
    }

    #[test]
    fn test_argmax_empty() {
        let data: Vec<f64> = Vec::new();
        assert_eq!(argmax(&data), None);
    }

    #[test]
    fn test_argmax_negatives() {
        let data = vec![-5.0, -4.0, -3.0, -2.0, -1.0];
        assert_eq!(argmax(&data), Some(4));
    }

    #[test]
    fn test_argmax_all_same() {
        let data = vec![2.0, 2.0, 2.0, 2.0];
        assert_eq!(argmax(&data), Some(0)); // It should return the first occurrence
    }
}
