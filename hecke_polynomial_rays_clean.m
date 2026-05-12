// ============================================================================
// hecke_polynomial_rays_clean.m
//
// Clean polynomial-ray calculator for the paper
//   "Growth in affine Hecke categories".
//
// This script computes nu_delta(b_w), i.e. the standard-basis coefficient sum
// at v=1, along the ray families used in the paper's polynomial-growth sections.
//
// Implemented blocks:
//   (A) Affine A~1:
//       the two alternating rays of length l. The paper proves
//           p_l = nu_delta(b_w) = 2l+1.
//       This block prints the computed values and checks that identity.
//
//   (B) Affine A~2 wall rays:
//       x_j = 123123... of length j,
//       y_j = x_{j-2} j,
//       z_j = x_{j-3} j,
//       with final labels read modulo 3.
//       The paper studies
//           p_j^x = nu_delta(b_{x_j}),
//           p_j^y = nu_delta(b_{y_j}),
//           p_j^z = nu_delta(b_{z_j}),
//       and proves that all three grow quadratically.
//
// Dependencies:
//   Joel Gibson's IHecke package, attached via
//       AttachSpec("IHecke/IHecke.spec");
// ============================================================================

SetColumns(0);
AttachSpec("IHecke/IHecke.spec");
IHeckeVersion();

// ----------------------------------------------------------------------------
// User parameters
// ----------------------------------------------------------------------------

RUN_A1 := true;
RUN_A2_WALL := true;

A1_MAX_LENGTH := 20;
A2_MAX_J := 40;

PRINT_A1_ROWS := true;
PRINT_A2_ROWS := true;
PRINT_FINAL_VECTORS := true;
PRINT_A2_QUADRATIC_RATIOS := true; // Prints p_j / j^2 for large-scale sanity.

Z := Integers();
R := RealField(20);

// ----------------------------------------------------------------------------
// Generic helper: nu_delta(b_w)
// ----------------------------------------------------------------------------

function CoeffAt1(poly)
    return Z!Evaluate(poly, 1);
end function;

function NuDeltaFromWord(word, H, C)
    if #word eq 0 then
        bw := C!1;
    else
        bw := C.word;
    end if;
    hw := H!bw;
    _, coeffs := Support(hw);
    if #coeffs eq 0 then
        return Z!0;
    end if;
    return &+[ Z | CoeffAt1(cc) : cc in coeffs ];
end function;

function SupportSizeStandardFromWord(word, H, C)
    if #word eq 0 then
        bw := C!1;
    else
        bw := C.word;
    end if;
    hw := H!bw;
    supp, _ := Support(hw);
    return #supp;
end function;

// ============================================================================
// Block A: affine A~1 alternating rays
// ============================================================================

function AlternatingWord(start, len)
    // start=1 gives 1,12,121,1212,...
    // start=2 gives 2,21,212,2121,...
    if (start ne 1) and (start ne 2) then
        error "AlternatingWord(start,len): start must be 1 or 2";
    end if;
    if len lt 0 then
        error "AlternatingWord(start,len): len must be nonnegative";
    end if;

    wd := [ Integers() | ];
    for k in [1..len] do
        if ((k-1) mod 2) eq 0 then
            Append(~wd, start);
        else
            Append(~wd, 3-start);
        end if;
    end for;
    return wd;
end function;

if RUN_A1 then
    W1<s1,s2> := CoxeterGroup(GrpFPCox, "A~1");
    HAlg1, H1, C1 := ShortcutIHeckeAlgebra(W1);
    LPoly1<v1> := BaseRing(HAlg1);

    printf "\n============================================================\n";
    printf "Affine A~1 polynomial rays\n";
    printf "Computing p_l = nu_delta(b_w) along the two alternating rays\n";
    printf "Expected identity: p_l = 2*l + 1\n";
    printf "A1_MAX_LENGTH = %o\n", A1_MAX_LENGTH;
    printf "============================================================\n\n";

    a1_ray_start1 := [ Z | ];
    a1_ray_start2 := [ Z | ];

    for ell in [0..A1_MAX_LENGTH] do
        w1 := AlternatingWord(1, ell);
        w2 := AlternatingWord(2, ell);

        p1 := NuDeltaFromWord(w1, H1, C1);
        p2 := NuDeltaFromWord(w2, H1, C1);
        expected := 2*ell + 1;

        Append(~a1_ray_start1, p1);
        Append(~a1_ray_start2, p2);

        if PRINT_A1_ROWS then
            printf "ell=%o  ray1=%o  p=%o  ray2=%o  p=%o  expected=%o",
                   ell, w1, p1, w2, p2, expected;
            if (p1 eq expected) and (p2 eq expected) then
                printf "  [OK]\n";
            else
                printf "  [WARNING]\n";
            end if;
        end if;
    end for;

    if PRINT_FINAL_VECTORS then
        printf "\nA~1 vector, start 1: %o\n", a1_ray_start1;
        printf "A~1 vector, start 2: %o\n", a1_ray_start2;
    end if;
end if;

// ============================================================================
// Block B: affine A~2 wall-ray families x_j, y_j, z_j
// ============================================================================

function Mod3Label(k)
    return ((k - 1) mod 3) + 1;
end function;

function XWord(j)
    if j lt 0 then
        error "XWord(j): need j >= 0";
    end if;
    return [ Mod3Label(k) : k in [1..j] ];
end function;

function YWord(j)
    // y_j = x_{j-2} j; meaningful for j >= 3.
    if j lt 3 then
        return [ Integers() | ];
    end if;
    return XWord(j-2) cat [ Mod3Label(j) ];
end function;

function ZWord(j)
    // z_j = x_{j-3} j; meaningful for j >= 4.
    if j lt 4 then
        return [ Integers() | ];
    end if;
    return XWord(j-3) cat [ Mod3Label(j) ];
end function;

if RUN_A2_WALL then
    W2<s1,s2,s3> := CoxeterGroup(GrpFPCox, "A~2");
    HAlg2, H2, C2 := ShortcutIHeckeAlgebra(W2);
    LPoly2<v2> := BaseRing(HAlg2);

    printf "\n============================================================\n";
    printf "Affine A~2 wall-ray polynomial data\n";
    printf "Families: x_j, y_j, z_j from the paper\n";
    printf "Computing p_j^family = nu_delta(b_word)\n";
    printf "A2_MAX_J = %o\n", A2_MAX_J;
    printf "============================================================\n\n";

    px := [ Z | ];
    py_indices := [ Integers() | ];
    py := [ Z | ];
    pz_indices := [ Integers() | ];
    pz := [ Z | ];

    for j in [1..A2_MAX_J] do
        xwd := XWord(j);
        pxj := NuDeltaFromWord(xwd, H2, C2);
        Append(~px, pxj);

        y_valid := j ge 3;
        z_valid := j ge 4;

        if y_valid then
            ywd := YWord(j);
            pyj := NuDeltaFromWord(ywd, H2, C2);
            Append(~py_indices, j);
            Append(~py, pyj);
        end if;

        if z_valid then
            zwd := ZWord(j);
            pzj := NuDeltaFromWord(zwd, H2, C2);
            Append(~pz_indices, j);
            Append(~pz, pzj);
        end if;

        if PRINT_A2_ROWS then
            printf "j=%o  x_j=%o  p_x=%o", j, xwd, pxj;
            if PRINT_A2_QUADRATIC_RATIOS then
                printf "  p_x/j^2=%o", (R!pxj)/((R!j)^2);
            end if;

            if y_valid then
                printf "  |  y_j=%o  p_y=%o", YWord(j), py[#py];
                if PRINT_A2_QUADRATIC_RATIOS then
                    printf "  p_y/j^2=%o", (R!py[#py])/((R!j)^2);
                end if;
            else
                printf "  |  y_j=(not used)";
            end if;

            if z_valid then
                printf "  |  z_j=%o  p_z=%o", ZWord(j), pz[#pz];
                if PRINT_A2_QUADRATIC_RATIOS then
                    printf "  p_z/j^2=%o", (R!pz[#pz])/((R!j)^2);
                end if;
            else
                printf "  |  z_j=(not used)";
            end if;

            printf "\n";
        end if;
    end for;

    if PRINT_FINAL_VECTORS then
        printf "\nA~2 x-vector, j=1..%o: %o\n", A2_MAX_J, px;
        printf "A~2 y-indices: %o\n", py_indices;
        printf "A~2 y-vector:  %o\n", py;
        printf "A~2 z-indices: %o\n", pz_indices;
        printf "A~2 z-vector:  %o\n", pz;
    end if;
end if;

quit;
