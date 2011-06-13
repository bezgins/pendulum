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

   
let x_diff a b d sigma delta t x_now = Array2.of_array float64 c_layout [|
    [| b *. x_now.{0,0} *. ((1. -. sigma) *. x_now.{2,0} -. delta *. x_now.{1,0})|];
    [| x_now.{0,0} *. (1. -.
                       (1. -. delta) *. x_now.{1,0} +.
                       sigma *. x_now.{2,0}) |];
    [| a *. (x_now.{1,0} -. d *. x_now.{0,0})|]
    |]
;;

let x_0 = Array2.of_array float64 c_layout [|
    [| 1. |];
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


let _ =
    let a = 7.
    and b = 0.4
    and d = 1.17 
    and sigma = 0.284
    and delta = 0.681 in
    let f_x' = x_diff a b d sigma delta in
    let result = mx_rk4_up f_x' x_0 0. 0.001 30. in
    List.map (print_result) result
    
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
