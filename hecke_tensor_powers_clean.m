// ============================================================================
// hecke_tensor_powers_clean.m
//
// Clean tensor-power calculator for the paper
//   "Growth in affine Hecke categories".
//
// Main task:
//   Fix an affine Hecke algebra and an element b, then compute
//
//       b^n = Sum_x a_x(n) b_x
//
//   in the Kazhdan-Lusztig basis, and print
//
//       nu_n        = Sum_x a_x(n) evaluated at v=1,
//       nu_delta_n  = the standard-basis coefficient sum at v=1,
//       support_n   = the number of KL basis terms with nonzero coeff@1.
//
// The paper's guiding example is obtained with the default parameters:
//       AFFINE_TYPE := "A~2";
//       WORD        := [1,2,3];
//       ELEMENT_MODE := "KL";
//
// Then b = b_{123}, and the output begins
//       1, 3, 17, 83, 472, ...
// for nu_n.
//
// Dependencies:
//   Joel Gibson's IHecke package, attached via
//       AttachSpec("IHecke/IHecke.spec");
//
// IMPORTANT DISTINCTION:
//   ELEMENT_MODE = "KL"     means b = C.WORD, i.e. the KL basis element b_w.
//   ELEMENT_MODE = "SIMPLE" means b = C.[i1]*...*C.[ik].
//   The paper's tensor-power sequence uses ELEMENT_MODE = "KL".
// ============================================================================

SetColumns(0);
AttachSpec("IHecke/IHecke.spec");
IHeckeVersion();

// ----------------------------------------------------------------------------
// User parameters
// ----------------------------------------------------------------------------

AFFINE_TYPE := "A~2";       // Supported here: "A~1" or "A~2".
WORD := [1,2,3];            // Reduced word used to define the element.
ELEMENT_MODE := "KL";       // "KL" or "SIMPLE".
NMAX := 12;                  // Number of tensor powers to compute.

PRINT_EACH_ROW := true;      // Print one summary line for every n.
PRINT_FINAL_VECTORS := true; // Print vectors at the end.
PRINT_KL_EXPANSIONS := false;// Useful only for small NMAX; expensive/noisy.

// Optional scalar-poly-exp diagnostic. This prints
//     nu_n * n^NORMALIZING_EXPONENT / beta^n.
// For the paper's default A~2 example b_{123}, use exponent 3/2.
PRINT_NORMALIZED_RATIO := true;
NORMALIZING_EXPONENT := 3/2;

// If PRINT_KL_EXPANSIONS is true and AFFINE_TYPE = "A~2", the code can label
// many summands in the paper's notation x_j, y_j, z_j, theta(m,n), etc.
USE_PAPER_LABELS_A2 := true;
LABEL_SEARCH_DEPTH := 3*NMAX + 10;

// ----------------------------------------------------------------------------
// Build the affine Hecke algebra
// ----------------------------------------------------------------------------

if AFFINE_TYPE eq "A~1" then
    W<s1,s2> := CoxeterGroup(GrpFPCox, "A~1");
    gens := [s1,s2];
elif AFFINE_TYPE eq "A~2" then
    W<s1,s2,s3> := CoxeterGroup(GrpFPCox, "A~2");
    gens := [s1,s2,s3];
else
    error "AFFINE_TYPE must be \"A~1\" or \"A~2\" in this clean script.";
end if;

CoxeterDiagram(W);
HAlg, H, C := ShortcutIHeckeAlgebra(W);
LPoly<v> := BaseRing(HAlg);
Z := Integers();
R := RealField(30);

// ----------------------------------------------------------------------------
// Generic coefficient statistics
// ----------------------------------------------------------------------------

function CoeffAt1(poly)
    return Z!Evaluate(poly, 1);
end function;

function SumKLCoeffsAt1(elt_in_C, C)
    eltC := C!elt_in_C;
    _, coeffs := Support(eltC);
    if #coeffs eq 0 then
        return Z!0;
    end if;
    return &+[ Z | CoeffAt1(cc) : cc in coeffs ];
end function;

function NumKLSummandsAt1(elt_in_C, C)
    eltC := C!elt_in_C;
    _, coeffs := Support(eltC);
    n := 0;
    for cc in coeffs do
        if CoeffAt1(cc) ne 0 then
            n +:= 1;
        end if;
    end for;
    return n;
end function;

function SumStandardCoeffsAt1(elt_in_C, H, C)
    eltH := H!(C!elt_in_C);
    _, coeffs := Support(eltH);
    if #coeffs eq 0 then
        return Z!0;
    end if;
    return &+[ Z | CoeffAt1(cc) : cc in coeffs ];
end function;

function SimpleKLProduct(word, C)
    elt := C!1;
    for i in word do
        elt *:= C.[i];
    end for;
    return C!elt;
end function;

// ----------------------------------------------------------------------------
// A~2 paper-label machinery, used only for optional readable expansions
// ----------------------------------------------------------------------------

function WordSeqOfBasisElt(W, y)
    try
        return Eltseq(W!(y));
    catch e1
        try
            return Eltseq(y);
        catch e2
            return [ Integers() | ];
        end try;
    end try;
end function;

function CoxEltFromWord(W, gens, word)
    x := Identity(W);
    for i in word do
        x *:= gens[i];
    end for;
    return x;
end function;

function Mod3Label(k)
    return ((k - 1) mod 3) + 1;
end function;

function WordToString(word)
    s := "";
    if #word eq 0 then
        return "e";
    end if;
    for i in [1..#word] do
        s cat:= Sprint(word[i]);
    end for;
    return s;
end function;

function RelabelWord(word, sigma)
    return [ sigma[word[i]] : i in [1..#word] ];
end function;

function SigmaTag(sigma)
    return Sprint(sigma[1]) cat Sprint(sigma[2]) cat Sprint(sigma[3]);
end function;

function DecorateLabel(label, sigma)
    if sigma eq [1,2,3] then
        return label;
    else
        return label cat "[" cat SigmaTag(sigma) cat "]";
    end if;
end function;

SigmasA2 := [
    [1,2,3],
    [1,3,2],
    [2,1,3],
    [2,3,1],
    [3,1,2],
    [3,2,1]
];

function XWord(j)
    if j lt 0 then
        error "XWord(j): need j >= 0";
    end if;
    return [ Mod3Label(k) : k in [1..j] ];
end function;

function YWord(j)
    // y_j = x_{j-2} j, with the last label read modulo 3.
    if j lt 3 then
        return [ Integers() | ];
    end if;
    return XWord(j-2) cat [ Mod3Label(j) ];
end function;

function ZWord(j)
    // z_j = x_{j-3} j, with the last label read modulo 3.
    if j lt 4 then
        return [ Integers() | ];
    end if;
    return XWord(j-3) cat [ Mod3Label(j) ];
end function;

function ThetaWord(m, n)
    // theta(m,n)=123...(2m+1)(2m+2)(2m+1)2m...(2m-2n+1),
    // with labels read modulo 3.
    if (m lt 0) or (n lt 0) or (n gt m) then
        error "ThetaWord(m,n): need m >= 0 and 0 <= n <= m";
    end if;

    word := [ Integers() | ];
    for k in [1..2*m+1] do
        Append(~word, Mod3Label(k));
    end for;
    Append(~word, Mod3Label(2*m+2));

    k := 2*m + 1;
    while k ge 2*m - 2*n + 1 do
        Append(~word, Mod3Label(k));
        k -:= 1;
    end while;

    return word;
end function;

function LabelWordA2(W, gens, word, maxSearch, sigmas)
    // Recognition is by equality in W, not literal word equality.
    target := CoxEltFromWord(W, gens, word);
    L := #word;

    if L eq 0 then
        return "e";
    end if;

    for sigma in sigmas do
        // Cyclic wall words x_j.
        model := RelabelWord(XWord(L), sigma);
        if target eq CoxEltFromWord(W, gens, model) then
            return DecorateLabel("x_" cat Sprint(L), sigma);
        end if;

        // Wall neighbors y_j and z_j.
        for j in [1..maxSearch] do
            model := RelabelWord(YWord(j), sigma);
            if (#model gt 0) and (target eq CoxEltFromWord(W, gens, model)) then
                return DecorateLabel("y_" cat Sprint(j), sigma);
            end if;

            model := RelabelWord(ZWord(j), sigma);
            if (#model gt 0) and (target eq CoxEltFromWord(W, gens, model)) then
                return DecorateLabel("z_" cat Sprint(j), sigma);
            end if;
        end for;

        // Beyond-the-wall words theta(m,n).
        for m in [0..maxSearch] do
            for n in [0..m] do
                model := RelabelWord(ThetaWord(m,n), sigma);
                if target eq CoxEltFromWord(W, gens, model) then
                    return DecorateLabel(
                        "theta(" cat Sprint(m) cat "," cat Sprint(n) cat ")",
                        sigma
                    );
                end if;
            end for;
        end for;

        // theta(m,n)i.
        for m in [0..maxSearch] do
            for n in [0..m] do
                for i in [1..3] do
                    model := RelabelWord(ThetaWord(m,n) cat [i], sigma);
                    if target eq CoxEltFromWord(W, gens, model) then
                        return DecorateLabel(
                            "theta(" cat Sprint(m) cat "," cat Sprint(n) cat ")" cat Sprint(i),
                            sigma
                        );
                    end if;
                end for;
            end for;
        end for;

        // i theta(m,n).
        for m in [0..maxSearch] do
            for n in [0..m] do
                for i in [1..3] do
                    model := RelabelWord([i] cat ThetaWord(m,n), sigma);
                    if target eq CoxEltFromWord(W, gens, model) then
                        return DecorateLabel(
                            Sprint(i) cat "theta(" cat Sprint(m) cat "," cat Sprint(n) cat ")",
                            sigma
                        );
                    end if;
                end for;
            end for;
        end for;

        // i theta(m,n) j.
        for m in [0..maxSearch] do
            for n in [0..m] do
                for i in [1..3] do
                    for j in [1..3] do
                        model := RelabelWord([i] cat ThetaWord(m,n) cat [j], sigma);
                        if target eq CoxEltFromWord(W, gens, model) then
                            return DecorateLabel(
                                Sprint(i) cat "theta(" cat Sprint(m) cat "," cat Sprint(n) cat ")" cat Sprint(j),
                                sigma
                            );
                        end if;
                    end for;
                end for;
            end for;
        end for;
    end for;

    return WordToString(word);
end function;

function LabelBasisEltA2(W, gens, y, maxSearch, sigmas)
    return LabelWordA2(W, gens, WordSeqOfBasisElt(W, y), maxSearch, sigmas);
end function;

procedure PrintKLExpansionAt1(elt_in_C, W, gens, C :
        UsePaperLabelsA2 := false,
        MaxSearch := 20,
        Sigmas := SigmasA2)

    eltC := C!elt_in_C;
    supp, coeffs := Support(eltC);

    if #supp eq 0 then
        printf "0\n";
        return;
    end if;

    first := true;
    for i in [1..#supp] do
        c := CoeffAt1(coeffs[i]);
        if c ne 0 then
            if UsePaperLabelsA2 then
                label := LabelBasisEltA2(W, gens, supp[i], MaxSearch, Sigmas);
            else
                label := WordToString(WordSeqOfBasisElt(W, supp[i]));
            end if;

            if first then
                first := false;
            else
                printf " + ";
            end if;

            if c eq 1 then
                printf "b_{%o}", label;
            else
                printf "%o*b_{%o}", c, label;
            end if;
        end if;
    end for;
    printf "\n";
end procedure;

// ----------------------------------------------------------------------------
// Choose the element b
// ----------------------------------------------------------------------------

if ELEMENT_MODE eq "KL" then
    b := C.WORD;
elif ELEMENT_MODE eq "SIMPLE" then
    b := SimpleKLProduct(WORD, C);
else
    error "ELEMENT_MODE must be \"KL\" or \"SIMPLE\".";
end if;

beta := SumStandardCoeffsAt1(b, H, C);

printf "\n============================================================\n";
printf "Tensor-power computation in affine type %o\n", AFFINE_TYPE;
printf "WORD = %o\n", WORD;
printf "ELEMENT_MODE = %o\n", ELEMENT_MODE;
printf "beta = nu_delta(b) = %o\n", beta;
printf "NMAX = %o\n", NMAX;
printf "============================================================\n\n";

nu_seq := [ Z | ];
nu_delta_seq := [ Z | ];
support_seq := [ Integers() | ];

pow := C!1;
for n in [1..NMAX] do
    pow := C!(pow * b);

    nu_n := SumKLCoeffsAt1(pow, C);
    nu_delta_n := SumStandardCoeffsAt1(pow, H, C);
    support_n := NumKLSummandsAt1(pow, C);

    Append(~nu_seq, nu_n);
    Append(~nu_delta_seq, nu_delta_n);
    Append(~support_seq, support_n);

    if PRINT_EACH_ROW then
        printf "n=%o  nu=%o  support=%o  nu_delta=%o  beta^n=%o",
               n, nu_n, support_n, nu_delta_n, beta^n;
        if nu_delta_n eq beta^n then
            printf "  [OK]";
        else
            printf "  [WARNING: nu_delta != beta^n]";
        end if;
        if PRINT_NORMALIZED_RATIO then
            ratio := (R!nu_n) * ((R!n)^(R!NORMALIZING_EXPONENT)) / (R!(beta^n));
            printf "  normalized=%o", ratio;
        end if;
        printf "\n";
    end if;

    if PRINT_KL_EXPANSIONS then
        printf "KL expansion at v=1 for n=%o:\n", n;
        PrintKLExpansionAt1(
            pow, W, gens, C :
            UsePaperLabelsA2 := (USE_PAPER_LABELS_A2 and AFFINE_TYPE eq "A~2"),
            MaxSearch := LABEL_SEARCH_DEPTH,
            Sigmas := SigmasA2
        );
        printf "\n";
    end if;
end for;

if PRINT_FINAL_VECTORS then
    printf "\nnu sequence:        %o\n", nu_seq;
    printf "support sequence:   %o\n", support_seq;
    printf "nu_delta sequence:  %o\n", nu_delta_seq;
end if;

quit;
