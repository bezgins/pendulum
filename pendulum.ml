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

let f_magnit b sigma delta x =
    let x' = x.{0,0}
    and y' = x.{1,0}
    and z' = x.{2,0} in
    Array2.of_array float64 c_layout [|
    [| (b *. x') *. ((1. -. sigma) *. z' -. (delta *. y'))|];
    [| x' *. (((delta -. 1.) *. y') +. (sigma *. z')) |];
    [| 0. |]|]
;;


let x_magnit_diff a b f u t x = 
    let linear = a @*. x
    and non_linear = f x
    and control = u &*.b in
    (linear @+. non_linear) @+. control
;;

let x_0 = Array2.of_array float64 c_layout [|
    [| 3. +. (atan (-5.)) |];
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
    [| 1.2 ; 0. ; 0. |];
    [| 0. ; 0. ; 0. |];
    [| 0. ; 0. ; 0. |] |]
;;

let p_r = Array2.of_array float64 c_layout [|
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

let print_u_x b = 
    match b with
    (x , z, y) ->
        Printf.printf "%f\t" x;
        for i = 0 to Array2.dim1 y - 1 do
            for j = 0 to Array2.dim2 y - 1 do
                Printf.printf " %f\t" y.{i,j}
           done;
        done;
        Printf.printf "%f" z;
        print_newline ();
;;


let u_func b r h c k_t = 
    let (t, k) = k_t in
(*    and b' = transpose b in
    let u'_0 = (1./.r) &*. b' *)
    let u'_1_0 = k @*. c in
    let u'_1 = h @-. u'_1_0 in
    let u' = (1./.r) *. (b ><. u'_1) in
    (t, u')
;;

let u_func_n b r h x k = 
    (*let b' = transpose b in
    let u' = (1./.r) &*. (b' @*. (h @-. (k @*. x))) in*)
    let u'_0 = k @*. x in
    let u'_1 = h @-. u'_0 in
    let u' = (1./.r) *. (b ><. u'_1) in
    u'
;;

let h_diff a b q r f z_f k x t h = 
    let b' = transpose b in
    let h_diff_0 = (1. /. r) &*. ((b @*. b') @*. k) in
    let h_diff_part_0 = a @-. h_diff_0 in 
    let h_diff_part_0' = (transpose h_diff_part_0) @*. h
    and h_diff_part_1 = q @*. (z_f t)
    and h_diff_part_2 = k @*. (f x) in
    (h_diff_part_2 @-. h_diff_part_0') @-. h_diff_part_1
;;

let x_new h diff_funct prev u = 
    let (_, _, x)::_ = prev
    and (t, u') = u in
    let x_new = mx_rk4_step (diff_funct u') x t h in
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
    let x_new = mx_rk4_step (diff_func u) x t step in
    let u_new = u_f h_n' x_new k_n' in
    (t +. step, u_new, x_new)::prev
;;

let pi = 4. *. (atan 1.);;

let z t = Array2.of_array float64 c_layout [|
    [| 3. +. (atan (t *. 10. -. 5.)) |];
    [| 0. |];
    [| 0. |] |]
;;

let x_next calc_step h_diff_f h_new_f u_next_f 
    x_next_f x_diff_f z_f q b r
    krev k k0 k1 x0 x_prev = 
    (* h_new = h_n(0)..h_n(1) *)
    let _::h_new = 
        List.fold_left2 (h_new_f (-.calc_step) h_diff_f)
                ((1., q @*. (z_f 1.))::[])
                ((1., k1)::krev) x_prev in 
    let (_, h_new_0)::_ = h_new in
    let u_new_0 = u_next_f b r h_new_0 x0 k0 in
    (* x_n = x_n(1) .. x_n(0) *)
    let x_n =
        List.fold_left2 (x_next_f calc_step (u_next_f b r) x_diff_f)
                ((0., u_new_0, x0)::[])
                h_new k in
    x_n
;;
(* TODO: change 2 rev_map2's and fold_left to fold_left2 *)
let metric_func a b = 
    let diff x y = 
        let (_,_,x') = x
        and (_,_,y') = y in
        x' @-. y'
    in
    let diffs = List.rev_map2 diff a b in
    let dots = List.rev_map2 (><.) diffs diffs in
    List.fold_left (+.) 0. dots
;;

let fx eps calc_step metric_f x_func x_0 = 
    let rec fx' i metric_prev x_p =
        let _::x = x_func x_p in
        let metric  = (metric_f x_p x) in
        match (metric_prev -. metric) < eps with
          true  -> Printf.printf "%d\n" i; x
        | false -> fx' (i+1) metric x
    in
    let _::x_1 = x_func x_0 in
    let metric_0 = (metric_f x_0 x_1) in
    fx' 0 metric_0 x_1
;;

let _ =
    let in_file = Pervasives.open_in "in.txt" in
    let step = Scanf.fscanf in_file "%f " (fun x -> x) 
    and eps  = Scanf.fscanf in_file "%f " (fun x -> x) 
    and a = Scanf.fscanf in_file "%f " (fun x -> x) 
    and b = Scanf.fscanf in_file "%f " (fun x -> x)
    and d = Scanf.fscanf in_file "%f " (fun x -> x) 
    and r_r = Scanf.fscanf in_file "%f " (fun x -> x)
    and sigma = Scanf.fscanf in_file "%f " (fun x -> x)
    and delta = Scanf.fscanf in_file "%f " (fun x -> x) 
    and _ = q_r.{0,0} <- Scanf.fscanf in_file "%f " (fun x -> x)
    and _ = p_r.{0,0} <- Scanf.fscanf in_file "%f " (fun x -> x) in
    let k_1 = p_r in
    let a_r = a_magnit a d
    and b_r = b_magnit
    and f_r = f_magnit b sigma delta in
    let x_magnit_d = x_magnit_diff a_r b_r f_r 
    and f_x'= x_magnit_diff a_r b_r f_r 0.
    and h_magnit_d = h_diff a_r b_r q_r r_r f_r z in
    (* k_r = k(0) :: ... :: k(1)  *)
    let k_r = mx_rk4_down (k_magnit_diff a_r b_r q_r r_r) k_1 1. step 0. in
    (* k_rev = k(1) :: .. :: k(0) *)
    let _::k_rev = List.rev k_r in
    let (_,k_0)::_ = k_r in
    let h_0 = p_r @*. z(1.) in
    (* u_0 = u_0(0)...u_0(1) *)
    let u_0 = List.rev_map (u_func b_r r_r h_0 x_0) k_rev in
    let (_, u0)::_ = u_0 in
    let x_0'= List.fold_left (x_new step x_magnit_d) ((0. , u0, x_0)::[]) u_0 in

    let x_21' = fx eps step metric_func 
                (x_next step h_magnit_d h_new u_func_n 
                x_new_n x_magnit_d z q_r b_r r_r 
                k_rev k_r k_0 k_1 x_0) x_0' 
    in

    List.map (print_u_x) x_21';

;;
