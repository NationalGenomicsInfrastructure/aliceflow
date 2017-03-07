#!/usr/bin/env nextflow
// join() test by Maxime

A = Channel.from(['1','2','3',file('a')], ['4','5','6',file('b')])

A = A.view {"A: $it"}
A = A.map {a, b, c, d -> d}

A = A.view {"map: $it"}
A = A.toList()

A = A.view {"toList: $it"}

process output {
    input: file(x) from A

    output: stdout result

    script:
    bams = x.collect{"-I $it"}.join(' ')
    """
    echo $bams
    """
}

result.subscribe { println "suscribe: $it" }
