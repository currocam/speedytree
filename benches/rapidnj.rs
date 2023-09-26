use birc_rapidnj::property_tests::random_additive_tree::{
    distance_matrix_from_tree, random_unrooted_binary_tree,
};
use birc_rapidnj::rapid_nj::rapid_nj;
use criterion::{black_box, criterion_group, criterion_main, Criterion};
fn criterion_benchmark(c: &mut Criterion) {
    let original_tree = random_unrooted_binary_tree(250);
    let d = distance_matrix_from_tree(original_tree.clone());
    c.bench_function("rapidnj 250", |b| {
        b.iter(|| rapid_nj(black_box(d.clone())).unwrap())
    });
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);
