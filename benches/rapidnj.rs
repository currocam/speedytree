use birc_rapidnj::property_tests::random_additive_tree::{
    distance_matrix_from_tree, random_unrooted_binary_tree,
};
use birc_rapidnj::algo::neighbor_joining;
use criterion::{black_box, criterion_group, criterion_main, Criterion};
fn criterion_benchmark(c: &mut Criterion) {
    let original_tree = random_unrooted_binary_tree(250);
    let d = distance_matrix_from_tree(original_tree.clone());
    c.bench_function("rapidnj 250", |b| {
        b.iter(|| neighbor_joining(black_box(d.clone()), 200).unwrap())
    });
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);
