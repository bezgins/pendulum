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

let k_magnit_diff a b q r t k_now = 
    let a' = transpose a
    and b' = transpose b in
    let k_diff1 = -1. &*. (k_now @*. a) 
    and k_diff2 = a' @*. k_now 
    and k_diff3 = (1./.r) &*. (((k_now @*. b) @*. b') @*. k_now) in
    ((k_diff1 @-. k_diff2) @+. k_diff3) @-. q
;;


let a_magnit a d = Array2.of_array float64 c_layout [|
    [| 0. ; 0. ; 0. |];
    [| 1. ; 0. ; 0. |];
    [| -.a*.d; a; 0. |] |]
;;

let b_magnit = Array2.of_array float64 c_layout [|
    [| 1. |];
    [| 0. |];
    [| 0. |]|]
;;

let f_magnit b sigma delta x u = Array2.of_array float64 c_layout [|
    [| b *. x.{0,0} *. ((1. -. sigma) *. x.{2,0} -. delta *. x.{1,0})|];
    [| x.{0,0} *. ((delta -. 1.) *. x.{1,0} +. sigma *. x.{2,0}) |];
    [| 0. |]|]
;;


let x_magnit_diff a b f u t x = 
    ((a @*. x) @+. (u &*. b)) @+. (f x u)
;;

let x_0 = Array2.of_array float64 c_layout [|
    [| atan (-2.) |];
    [| 1. |];
    [| 1. |]|]
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

let solve' k g k0 g0 t0 h t_end =
    let rec rk4_rec result k g k_now g_now t h t_end =
        let k_new = mx_rk4_step k k_now t (h) 
        and g_new = mx_rk4_step (g k_now) g_now t (h)
        and t_new = t +. h in
        match t_new > t_end with
          true  -> ((t_new, k_new, g_new)::result)
        | false -> rk4_rec ((t_new, k_new, g_new)::result) k g k_new g_new t_new h t_end
    in
    rk4_rec [] k g k0 g0 t0 h t_end
;;

let q_r = Array2.of_array float64 c_layout [|
    [| 1. ; 0. ; 0. |];
    [| 0. ; 0. ; 0. |];
    [| 0. ; 0. ; 0. |] |]
;;


let print_result b = 
    match b with
    (x , y) ->
        Printf.printf "%f\t" x;
        for i = 0 to Array2.dim1 y - 1 do
            for j = 0 to Array2.dim2 y - 1 do
                Printf.printf " %f\t" y.{i,j}
           done;
        done;
(*        Printf.printf "%f" z;*)
        print_newline ();
;;

let u_func b r h x k_t = 
    let (t, k) = k_t
    and b' = transpose b in
    let u' = (1./.r) &*. (b' @*. (h @-. (k @*. x))) in
    (t, u')
;;

let h_diff a b k q r f x z h = 
    let b' = transpose b in
    let h_diff_0 = (1. /. r) &*. ((b @*. b') @*. k) in
    let h_diff_part_0 = a @-. h_diff_0 in 
    let h_diff_part_0' = transpose h_diff_part_0
    and h_diff_part_1 = q @*. z
    and h_diff_part_2 = k @*. (f x) in
    h_diff_part_2 @-. h_diff_part_0' @-. h_diff_part_1
;;

let z t = Array2.of_array float64 c_layout [|
    [| atan ((t *. 10. -. 5.) +. 3.) |];
    [| 0. |];
    [| 0. |] |]
;;


let _ =
    let a = 7.
    and b = 0.4
    and d = 1.17 
    and r_r = 0.5
    and k_0 = q_r
    and sigma = 0.284
    and delta = 0.681 in
    let a_r = a_magnit a d
    and b_r = b_magnit
    and f_r = f_magnit b sigma delta in
    let f_x' = x_magnit_diff a_r b_r f_r 0. in
    let k_r = mx_rk4_down (k_magnit_diff a_r b_r q_r r_r) k_0 1. 0.001 0. in
    let h_0 = q_r @*. z(1.) in
    let u_0 = List.map (u_func b_r r_r h_0 x_0) k_r in
    List.map (print_result) u_0
    (*
    let result = mx_rk4_up f_x' x_0 0. 0.001 30. in
    List.map (print_result) result
    *)
    
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
