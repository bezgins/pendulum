(*
Copyright (c) 2011, Svyatoslav Bezgin <bezgin@rcbd.org> All rights reserved.

Redistribution and use in source and binary forms, with or without modification,
are permitted provided that the following conditions are met:

 *  Redistributions of source code must retain the above copyright notice,
    this list of conditions and the following disclaimer.
 *  Redistributions in binary form must reproduce the above copyright notice,
    this list of conditions and the following disclaimer in the documentation
    and/or other materials provided with the distribution.
 *  Neither the name of the Tambov State Technical University nor the names of
    its contributors may be used to endorse or promote products derived from this
    software without specific prior written permission.
            
THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY
EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT
SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT
OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*)

open Batteries;;
open Enum;;
open Print;;
open Bigarray;;
open Matrix;;

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
