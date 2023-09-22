use std::io::{self, Read, Write};
use std::process::{Command, Stdio};

#[test]
fn wikipedia_example_rapid() {
    //https://en.wikipedia.org/wiki/Neighbor_joining
    let input = "5
    a	0	5	9	9	8
    b	5	0	10	10	9
    c	9	10	0	8	7
    d	9	10	8	0	3
    e	8	9	7	3	0
";

    let expected_output = "((c:4.0,(d:2.0,e:1.0):2.0):3.0,a:2.0,b:3.0);";

    let mut child = Command::new("target/debug/birc-rapidnj")
        .arg("--rapid")
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .spawn()
        .expect("Failed to spawn child process");

    let stdin = child.stdin.as_mut().unwrap();
    stdin.write_all(input.as_bytes()).unwrap();

    let mut output = String::new();
    child
        .stdout
        .as_mut()
        .unwrap()
        .read_to_string(&mut output)
        .unwrap();

    let status = child.wait().unwrap();
    assert!(status.success());

    assert_eq!(output.trim(), expected_output);
}

#[test]
fn simple_tree_rapidnj() {
    let input = "6
    Mouse     0.0000 1.5232 1.4841 1.4465 1.4389 1.4629 
    Gibbon    1.5232 0.0000 0.7115 0.5958 0.6179 0.5583 
    Orang     1.4841 0.7115 0.0000 0.4631 0.5061 0.4710 
    Gorilla   1.4465 0.5958 0.4631 0.0000 0.3484 0.3083 
    Chimp     1.4389 0.6179 0.5061 0.3484 0.0000 0.2692 
    Human     1.4629 0.5583 0.4710 0.3083 0.2692 0.0000
";
    let expected_output = "(((Gorilla:0.158225,(Chimp:0.15009999999999988,Human:0.11910000000000012):0.03552500000000003):0.03500000000000006,Orang:0.27664999999999997):0.05954999999999988,Mouse:1.1802124999999999,Gibbon:0.3429875000000002);";
    // Parse the output so only one decimal place is shown

    let mut child = Command::new("target/debug/birc-rapidnj")
        .arg("--rapid")
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .spawn()
        .expect("Failed to spawn child process");

    let stdin = child.stdin.as_mut().unwrap();
    stdin.write_all(input.as_bytes()).unwrap();

    let mut output = String::new();
    child
        .stdout
        .as_mut()
        .unwrap()
        .read_to_string(&mut output)
        .unwrap();

    let status = child.wait().unwrap();
    assert!(status.success());

    assert_eq!(output.trim(), expected_output);
}
