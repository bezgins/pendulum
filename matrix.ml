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

let mx_rk4_step f y t h = 
    let step1 = f t y in
    let step2 = f (t +. (h /. 2.)) (y @+. ((h /. 2.) &*. step1)) in
    let step3 = f (t +. (h /. 2.)) (y @+. ((h /. 2.) &*. step2)) in
    let step4 = f (t +. h) (y @+. (h &*. step3)) in
    y @+. ((h /. 6.) &*.
        (((step1 @+. (2.0 &*. step2)) @+. (2.0 &*. step3)) @+. step4 ))
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

