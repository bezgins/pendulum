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

let f_magnit b sigma delta x = Array2.of_array float64 c_layout [|
    [| b *. x.{0,0} *. ((1. -. sigma) *. x.{2,0} -. delta *. x.{1,0})|];
    [| x.{0,0} *. ((delta -. 1.) *. x.{1,0} +. sigma *. x.{2,0}) |];
    [| 0. |]|]
;;


let x_magnit_diff a b f u t x = 
    ((a @*. x) @+. (u &*. b)) @+. (f x)
;;

let x_0 = Array2.of_array float64 c_layout [|
    [| 3. (*(atan (-5.)) +. 3.*) |];
    [| 4. |];
    [| 4. |]|]
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
    [| 10. ; 0. ; 0. |];
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

let print_u_x b = 
    match b with
    (x , z, y) ->
        Printf.printf "%f\t" x;
        for i = 0 to Array2.dim1 y - 1 do
            for j = 0 to Array2.dim2 y - 1 do
                Printf.printf " %f\t" y.{i,j}
           done;
        done;
        Printf.printf "%f" z.{0,0};
        print_newline ();
;;


let u_func b r h x k_t = 
    let (t, k) = k_t
    and b' = transpose b in
    let u' = (1./.r) &*. (b' @*. (h @-. (k @*. x))) in
    (t, u')
;;

let u_func_n b r h x k = 
    let b' = transpose b in
    let u' = (1./.r) &*. (b' @*. (h @-. (k @*. x))) in
    u'
;;

let h_diff a b q r f z_f k x t h = 
    let b' = transpose b in
    let h_diff_0 = (1. /. r) &*. ((b @*. b') @*. k) in
    let h_diff_part_0 = a @-. h_diff_0 in 
    let h_diff_part_0' = transpose h_diff_part_0
    and h_diff_part_1 = q @*. (z_f t)
    and h_diff_part_2 = k @*. (f x) in
    (h_diff_part_2 @-. h_diff_part_0') @-. h_diff_part_1
;;

let x_new h diff_funct prev u = 
    let (t, _, x)::_ = prev
    and (_, u') = u in
    let x_new = mx_rk4_step (diff_funct u'.{0,0}) x t h in
    (t +. h, u', x_new)::prev
;;

let h_new step diff_funct prev k_r x_p = 
    let (t, k) = k_r
    and (_, _, x) = x_p
    and (_, h)::_ = prev in
    let h_new =  mx_rk4_step (diff_funct k x) h t step in
    (t +. step, h_new)::prev
;;

let x_new_n step u_f diff_func prev h_n k_n =
    let (t, u, x)::_  = prev 
    and (_, h_n') = h_n
    and (_, k_n') = k_n in
    let x_new = mx_rk4_step (diff_func u.{0,0}) x t step in
    let u_new = u_f h_n' x_new k_n' in
    (t +. step, u_new, x_new)::prev
;;

let z t = Array2.of_array float64 c_layout [|
    [| 3. +. t(*(atan (t *. 10. -. 5.)) +. 3.*) |];
    [| 0. |];
    [| 0. |] |]
;;

let x_next h_diff_f h_new_f u_next_f 
    x_next_f x_diff_f z_f q b r
    krev k k0 k1 x0 x_prev = 
    let _::h_new = 
        List.fold_left2 (h_new_f (-.0.001) h_diff_f)
                ((1., q @*. (z_f 1.))::[])
                ((1., k1)::krev) x_prev in 
    let (_, h_new_0)::_ = h_new in
    let u_new_0 = u_next_f b r h_new_0 x0 k0 in
    let x_n =
        List.fold_left2 (x_next_f 0.001 (u_next_f b r) x_diff_f)
                ((0., u_new_0, x0)::[])
                h_new k in
    x_n
;;

let _ =
    let a = 7.
    and b = 0.4
    and d = 1.17 
    and r_r = 1.
    and k_1 = q_r
    and sigma = 0.284
    and delta = 0.681 in
    let a_r = a_magnit a d
    and b_r = b_magnit
    and f_r = f_magnit b sigma delta in
    let x_magnit_d = x_magnit_diff a_r b_r f_r 
    and f_x'= x_magnit_diff a_r b_r f_r 0.
    and h_magnit_d = h_diff a_r b_r q_r r_r f_r z in
    let k_r = mx_rk4_down (k_magnit_diff a_r b_r q_r r_r) k_1 1. 0.001 0. in
    let _::k_rev = List.rev k_r 
    and (_,k_0)::_ = k_r in
    let h_0 = q_r @*. z(1.) in
    let u_0 = List.rev_map (u_func b_r r_r h_0 x_0) k_rev in
    let (_, u0)::_ = u_0 in
    let x_0'= List.fold_left (x_new 0.001 x_magnit_d) ((0. , u0, x_0)::[]) u_0 in

    let _::x_1 = x_next h_magnit_d h_new u_func_n 
                x_new_n x_magnit_d z q_r b_r r_r 
                k_rev k_r k_0 k_1 x_0 x_0' in

    let _::x_2 = x_next h_magnit_d h_new u_func_n 
                x_new_n x_magnit_d z q_r b_r r_r 
                k_rev k_r k_0 k_1 x_0 x_1 in

    let _::x_3 = x_next h_magnit_d h_new u_func_n 
                x_new_n x_magnit_d z q_r b_r r_r 
                k_rev k_r k_0 k_1 x_0 x_2 in

    let _::x_4 = x_next h_magnit_d h_new u_func_n 
                x_new_n x_magnit_d z q_r b_r r_r 
                k_rev k_r k_0 k_1 x_0 x_3 in

    let _::x_5 = x_next h_magnit_d h_new u_func_n 
                x_new_n x_magnit_d z q_r b_r r_r 
                k_rev k_r k_0 k_1 x_0 x_4 in

    let _::x_6 = x_next h_magnit_d h_new u_func_n 
                x_new_n x_magnit_d z q_r b_r r_r 
                k_rev k_r k_0 k_1 x_0 x_5 in
    let _::x_7 = x_next h_magnit_d h_new u_func_n 
                x_new_n x_magnit_d z q_r b_r r_r 
                k_rev k_r k_0 k_1 x_0 x_6 in

    let _::x_8 = x_next h_magnit_d h_new u_func_n 
                x_new_n x_magnit_d z q_r b_r r_r 
                k_rev k_r k_0 k_1 x_0 x_7 in

    let _::x_9 = x_next h_magnit_d h_new u_func_n 
                x_new_n x_magnit_d z q_r b_r r_r 
                k_rev k_r k_0 k_1 x_0 x_8 in

    let _::x_10 = x_next h_magnit_d h_new u_func_n 
                x_new_n x_magnit_d z q_r b_r r_r 
                k_rev k_r k_0 k_1 x_0 x_9 in

    let _::x_11 = x_next h_magnit_d h_new u_func_n 
                x_new_n x_magnit_d z q_r b_r r_r 
                k_rev k_r k_0 k_1 x_0 x_10 in

    let _::x_12 = x_next h_magnit_d h_new u_func_n 
                x_new_n x_magnit_d z q_r b_r r_r 
                k_rev k_r k_0 k_1 x_0 x_11 in

    let _::x_13 = x_next h_magnit_d h_new u_func_n 
                x_new_n x_magnit_d z q_r b_r r_r 
                k_rev k_r k_0 k_1 x_0 x_12 in

    let _::x_14 = x_next h_magnit_d h_new u_func_n 
                x_new_n x_magnit_d z q_r b_r r_r 
                k_rev k_r k_0 k_1 x_0 x_13 in

    let _::x_15 = x_next h_magnit_d h_new u_func_n 
                x_new_n x_magnit_d z q_r b_r r_r 
                k_rev k_r k_0 k_1 x_0 x_14 in

    let _::x_16 = x_next h_magnit_d h_new u_func_n 
                x_new_n x_magnit_d z q_r b_r r_r 
                k_rev k_r k_0 k_1 x_0 x_15 in
    let _::x_17 = x_next h_magnit_d h_new u_func_n 
                x_new_n x_magnit_d z q_r b_r r_r 
                k_rev k_r k_0 k_1 x_0 x_16 in

    let _::x_18 = x_next h_magnit_d h_new u_func_n 
                x_new_n x_magnit_d z q_r b_r r_r 
                k_rev k_r k_0 k_1 x_0 x_17 in

    let _::x_19 = x_next h_magnit_d h_new u_func_n 
                x_new_n x_magnit_d z q_r b_r r_r 
                k_rev k_r k_0 k_1 x_0 x_18 in

    let _::x_20 = x_next h_magnit_d h_new u_func_n 
                x_new_n x_magnit_d z q_r b_r r_r 
                k_rev k_r k_0 k_1 x_0 x_19 in

    let _::x_21 = x_next h_magnit_d h_new u_func_n 
                x_new_n x_magnit_d z q_r b_r r_r 
                k_rev k_r k_0 k_1 x_0 x_20 in


    List.map (print_u_x) x_0';
    List.map (print_u_x) x_10;
    List.map (print_u_x) x_21;
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
