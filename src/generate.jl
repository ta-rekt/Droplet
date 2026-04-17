# . - -- .- .-.- . .- .-- -.- .- .- --.-.... - -.-
# TODO: Finish this shit
export generate
function generate(n, Q, idist)

    s, t = zeros(Int, n), zeros(n)
    
    s[1] = sample(eachindex(idist), Weights(idist))
    q = -Q[s[1], s[1]]

    _DD = _ddist(Q)

    for i in 1:n-1
        si = s[i]
        q = -Q[si, si]

        # t
        e = Exponential(1/q)
        t[i] = rand(e)

        # s        
        _dd = _DD[si, :]
        s[i+1] = sample(eachindex(_dd), Weights(_dd))
    end

    # t
    si = s[end]
    q = -Q[si, si]
    e = Exponential(1/q)
    t[end] = rand(e)

    return s, t
end

export generate!
function generate!(m::LLROptModel)
    # TODO: inline
    s, t = generate(m.n, m.Q_c, m.idist)
    
    m.s = s
    m.t = t
    return m
end