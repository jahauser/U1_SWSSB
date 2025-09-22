SWAP = [1 0 0 0
        0 0 1 0
        0 1 0 0
        0 0 0 1]

PauliZ = [1 0
          0 -1]

ITensors.op(::OpName"a",::SiteType"Qubit") = [0 1
                                              0 0]
ITensors.op(::OpName"a†",::SiteType"Qubit") = [0 0
                                              1 0]