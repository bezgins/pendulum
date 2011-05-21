open Batteries;;
open Enum;;
open Print;;
open Bigarray;;

(* Matrix transpose *)
let transpose b = 
    let dim1    = Array2.dim1 b
    and dim2    = Array2.dim2 b in
    let kind    = Array2.kind b
    and layout  = Array2.layout b in
    let b' = Array2.create kind layout dim2 dim1 in
    for i = 0 to pred dim1 do
        for j = 0 to pred dim2 do
            b'.{j,i} <- b.{i,j}
        done;
    done;
    (b')
;;

(* Generic matrix multiplication *)
let mx_mult plus mul nul x y =
    let x0 = Array2.dim1 x
    and y0 = Array2.dim1 y in
    let y1 = if y0 = 0 then 0 else Array2.dim2 y in
    let kind = Array2.kind x
    and layout = Array2.layout x in
    let z = Array2.create kind layout x0 y1 in
    for i = 0 to x0 - 1 do
        for j = 0 to y1 - 1 do
            z.{i,j} <- nul;
            for k = 0 to y0 - 1 do
                z.{i,j} <- plus z.{i,j} (mul x.{i,k} y.{k,j})
            done
        done
    done;
    z
;;

(* Float matrix multiplication *)
let (@*.) a b = 
    mx_mult (+.) ( *. ) 0. a b
;;

(* Generic matrix negation *)
let mx_neg nul minus b = 
    let dim1    = Array2.dim1 b
    and dim2    = Array2.dim2 b in
    let kind    = Array2.kind b
    and layout  = Array2.layout b in
    let b' = Array2.create kind layout dim1 dim2 in
     for i = 0 to pred dim1 do
        for j = 0 to pred dim2 do
            b'.{i,j} <- minus nul b.{i,j}
        done;
    done;
    (b')
;;

(* float matrix negation *)
let (~@-.) a = 
    mx_neg 0. (-.) a
;;

(* Generic matrix addition *)
let mx_add plus a b = 
    let dim1    = Array2.dim1 a
    and dim2    = Array2.dim2 a in
    let kind    = Array2.kind a
    and layout  = Array2.layout a in
    let b' = Array2.create kind layout dim1 dim2 in
     for i = 0 to pred dim1 do
        for j = 0 to pred dim2 do
            b'.{i,j} <- plus a.{i,j} b.{i,j}
        done;
    done;
    (b')
;;

(* Float matrix addition *)
let (@+.) a b =
    mx_add (+.) a b
;;

(* Float matrix substraction *)
let (@-.) a b =
    mx_add (-.) a b
;;


(* Generic scalar multiplication *)
let mx_scalar_mul mul a b = 
    let dim1    = Array2.dim1 b
    and dim2    = Array2.dim2 b in
    let kind    = Array2.kind b
    and layout  = Array2.layout b in
    let b' = Array2.create kind layout dim1 dim2 in
     for i = 0 to pred dim1 do
        for j = 0 to pred dim2 do
            b'.{i,j} <- mul a b.{i,j}
        done;
    done;
    (b')
;;

(* Float scalar multiplication*)
let (&*.) a b =
    mx_scalar_mul ( *. ) a b
;;

let array2_display print newline b = 
    for i = 0 to Array2.dim1 b - 1 do
        for j = 0 to Array2.dim2 b - 1 do
            print b.{i,j}
        done;
        newline();
    done;
;;

let some_print b = 
    match b with
    (x , y, z) ->
        Printf.printf "%f\t" x;
        for i = 0 to Array2.dim1 y - 1 do
            for j = 0 to Array2.dim2 y - 1 do
                Printf.printf " %f\t" y.{i,j}
           done;
        done;
        for i = 0 to Array2.dim1 z - 1 do
            for j = 0 to Array2.dim2 z - 1 do
                Printf.printf " %f\t" z.{i,j}
           done;
        done;
        print_newline ();
;;

let rk4_step f y t h = 
    let k1 = f t y in
    let k2 = f (t +. 0.5 *. h) (y +. 0.5 *. h *. k1) in
    let k3 = f (t +. 0.5 *. h) (y +. 0.5 *. h *. k2) in
    let k4 = f (t +. h) (y +. h *. k3) in
    y +. h *. (k1 +. 2.0 *. k2 +. 2.0 *. k3 +. k4) /. 6.0
;;

let mx_rk4_step f y t h = 
    let step1 = f t y in
    let step2 = f (t +. (h /. 2.)) (y @+. ((h /. 2.) &*. step1)) in
    let step3 = f (t +. (h /. 2.)) (y @+. ((h /. 2.) &*. step2)) in
    let step4 = f (t +. h) (y @+. (h &*. step3)) in
    y @+. ((h /. 6.) &*.
        (((step1 @+. (2.0 &*. step2)) @+. (2.0 &*. step3)) @+. step4 )
        )
;;

let rk4 f y0 t0 h t_end = 
    let rec rk4_rec result f y t h t_end =
        let y_new = rk4_step f y t h in
        let t_new = t +. h in
        match t_new > t_end with
          true  -> (t_new, y_new) :: result
        | false -> rk4_rec ((t_new, y_new)::result) f y_new t_new h t_end
    in
    rk4_rec [] f y0 t0 h t_end
;;

let mx_rk4_up f y0 t0 h t_end =
    let rec rk4_rec result f y t h t_end =
        let y_new = mx_rk4_step f y t h in
        let t_new = t +. h in
        match t_new > t_end with
          true  -> (t_new, y_new) :: result
        | false -> rk4_rec ((t_new, y_new)::result) f y_new t_new h t_end
    in
    rk4_rec [] f y0 t0 h t_end
;;

let mx_rk4_down f y0 t0 h t_end =
    let rec rk4_rec result f y t h t_end =
        let y_new = mx_rk4_step f y t (-.h) in
        let t_new = t -. h in
        match t_new < t_end with
          true  -> result
        | false -> rk4_rec ((t_new, y_new)::result) f y_new t_new h t_end
    in
    rk4_rec [] f y0 t0 h t_end
;;

let pi = 4.0 *. atan 1.0;;

let a_r eps t = Array2.of_array float64 c_layout [|
    [| 0. ; 1. |];
    [| -.eps; -. cos (sin (2. *. pi *. t)) |] |]
;;

let b_r = Array2.of_array float64 c_layout [|
    [| 0. |];
    [| 1. |] |]
;;


let h_r t =  Array2.of_array float64 c_layout [|
    [| 0. |];
    [|(cos (sin (2. *. pi *. t))) *. (sin (2. *. pi *. t)) -.
       sin (sin (2. *. pi *. t))  |] |]
;;

let q_r = Array2.of_array float64 c_layout [|
    [| 5. ; 0.1|];
    [| 0.1 ; 10.|] |]
;;

let f_r = Array2.of_array float64 c_layout [|
    [| 1. ; 0.|];
    [| 0. ; 1.|] |]
;;

let phi_r t = Array2.of_array float64 c_layout [|
    [| sin (2. *. pi *. t) |] ;
    [|2. *. pi *. cos (2. *. pi *. t) |] |]
;;

let k_diff a b q r t k_now = 
    let at = a t in
    let at' = transpose at
    and b' = transpose b in
    let k_diff1 = -1. &*. (k_now @*. at) 
    and k_diff2 = at' @*. k_now 
    and k_diff3 = (1./.r) &*. (((k_now @*. b) @*. b') @*. k_now) in
    ((k_diff1 @-. k_diff2) @+. k_diff3) @-. q
;;

let g_diff a phi h b q r k t g_now = 
    let at = a t 
    and b' = transpose b 
    and phit = phi t
    and ht = h t in
    let g_diff_part = (b @*. b') @*. k in
    let g_diff1 = (transpose (at @-. g_diff_part)) @*. g_now
    and g_diff2 = q @*. phit
    and g_diff3 = k @*. ht in
    ((-.1./.r) &*. g_diff1) @-. g_diff2 @+. g_diff3
;;

let ustar b r g_now k_now x_now = 
    let b' = transpose b
    and u_star_diff_1 = k_now @*. x_now in
    let u_star_diff_2 = g_now @-. u_star_diff_1 in
    (1./.r) &*. (b' @*. u_star_diff_2)
;;
    
let x_diff eps u t x_now  = Array2.of_array float64 c_layout [|
    [| x_now.{1,0} |];
    [| -.sin(x_now.{0,0}) -. (eps *. x_now.{1,0}) +. u |] |]
;;

let x_r = Array2.of_array float64 c_layout [|
    [| 0. |];
    [| 0. |] |]
;;

let solve k g k0 g0 t0 h t_end =
    let rec rk4_rec result k g k_now g_now t h t_end =
        let k_new = mx_rk4_step k k_now t (-.h) 
        and g_new = mx_rk4_step (g k_now) g_now t (-.h)
        and t_new = t -. h in
        match t_new < t_end with
          true  -> ((t_new, k_new, g_new)::result)
        | false -> rk4_rec ((t_new, k_new, g_new)::result) k g k_new g_new t_new h t_end
    in
    rk4_rec [] k g k0 g0 t0 h t_end
;;

let solve_x u_func x_func k_g x0 t0 h t_end = 
    let rec solve' result k_g_now x_now t h t_end = 
        let (_, k_now, g_now)::tail = k_g_now in
        let u = u_func g_now k_now x_now in
        let x_new = mx_rk4_step (x_func u.{0,0}) x_now t h
        and t_new = t +. h in
        match t_new > t_end with
          true  -> (t_new, x_new, u.{0,0}) :: result
        | false -> solve' ((t_new, x_new, u.{0,0})::result) tail x_new t_new h t_end
    in
    solve' [] k_g x0 t0 h t_end
;;

let print_result b = 
    match b with
    (x , y, z) ->
        Printf.printf "%f\t" x;
        for i = 0 to Array2.dim1 y - 1 do
            for j = 0 to Array2.dim2 y - 1 do
                Printf.printf " %f\t" y.{i,j}
           done;
        done;
        Printf.printf "%f" z;
        print_newline ();
;;


let _ =
    let eps = 0.1
    and r = 0.0074184 in
    let k = (k_diff (a_r eps) b_r q_r r)
    and g = g_diff (a_r eps) phi_r h_r b_r q_r r
    and g_start = f_r @*. (phi_r 1.) in
    let k_g = solve k g f_r g_start 1.0 0.0001 0.0 in
    let result = solve_x (ustar b_r r) (x_diff eps) k_g x_r 0. 0.0001 1. in
    let (_, x1,_)::tail = result in
    let result2 = solve_x (ustar b_r r) (x_diff eps) k_g x1 1. 0.0001 1.9999 in
    let (_, x2,_)::tail = result2 in
    let result3 = solve_x (ustar b_r r) (x_diff eps) k_g x2 2. 0.0001 2.9999 in
    let (_, x3,_)::tail = result3 in
    let result4 = solve_x (ustar b_r r) (x_diff eps) k_g x3 3. 0.0001 3.9999 in
    let (_, x4,_)::tail = result4 in
    let result5 = solve_x (ustar b_r r) (x_diff eps) k_g x4 4. 0.0001 4.9999 in
    let (_, x5,_)::tail = result5 in
    let result6 = solve_x (ustar b_r r) (x_diff eps) k_g x5 5. 0.0001 5.9999 in
    let (_, x6,_)::tail = result6 in
    let result7 = solve_x (ustar b_r r) (x_diff eps) k_g x6 6. 0.0001 6.9999 in
    List.map (print_result) (List.rev result);
    List.map (print_result) (List.rev result2);
    List.map (print_result) (List.rev result3);
    List.map (print_result) (List.rev result4);
    List.map (print_result) (List.rev result5);
    List.map (print_result) (List.rev result6);
    List.map (print_result) (List.rev result7);

(*    List.map (some_print) result *)
(*    let k_diff = mx_rk4_step k f_r 1.0 (-0.0005) in
    let g_diff = mx_rk4_step (g f_r) g_start 1.0 (-0.0005) in
*)
    (*
    let res = mx_rk4_down k f_r 1.0 0.0005 0.0 in
    List.map (some_print) res
    array2_display (Printf.printf " %f") print_newline g_diff
    *)
;;
