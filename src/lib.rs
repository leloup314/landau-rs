mod pdfs;


#[cfg(test)]
mod tests {
    use super::pdfs::langau_pdf;
    #[test]
    fn it_works() {
        let result = 2 + 2;
        assert_eq!(result, 4);
        let a = 15f64;
        let b = 1f64;
        let c = 1f64;
        let d = 1f64;
        assert_eq!(langau_pdf(a, b, c, d), 0.006297828887560015);
    }
}
