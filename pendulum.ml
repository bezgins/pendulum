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

(* Число П. На всякий случай *)
let pi = 4. *. (atan 1.);;

(* Матричное уравнение Риккати *)
let k_magnit_diff a b q r t k_now = 
    let a' = transpose a
    and b' = transpose b in
    let k_diff1 = -1. &*. (k_now @*. a) 
    and k_diff2 = a' @*. k_now 
    and k_diff3 = (1./.r) &*. (((k_now @*. b) @*. b') @*. k_now) in
    ((k_diff1 @-. k_diff2) @+. k_diff3) @-. q
;;

(* Матрица A модели X'= AX + Bu + f(X) *)
let a_magnit a d = Array2.of_array float64 c_layout [|
    [| 0. ; 0. ; 0. |];
    [| 1. ; 0. ; 0. |];
    [| -.a*.d; a; 0. |] |]
;;

(* Матрица B модели *)
let b_magnit = Array2.of_array float64 c_layout [|
    [| 1. |];
    [| 0. |];
    [| 0. |]|]
;;

(* Матричная функция f модели *)
let f_magnit b sigma delta x =
    let x' = x.{0,0}
    and y' = x.{1,0}
    and z' = x.{2,0} in
    Array2.of_array float64 c_layout [|
    [| (b *. x') *. ((1. -. sigma) *. z' -. (delta *. y'))|];
    [| x' *. (((delta -. 1.) *. y') +. (sigma *. z')) |];
    [| 0. |]|]
;;

(* Собственно само дифф. уравнение модели Магницкого с управлением *)
let x_magnit_diff a b f u t x = 
    let linear = a @*. x
    and non_linear = f x
    and control = u &*.b 
    and p = 
        if  ((t > 0.199) && (t < 0.201)) ||
            ((t > 0.399) && (t < 0.401)) ||
            ((t > 0.601) && (t < 0.603)) ||
            ((t > 0.801) && (t < 0.803))
        then
            200.
        else
            if  ((t > 0.201) && (t < 0.203)) || 
                ((t > 0.401) && (t < 0.403)) ||
                ((t > 0.599) && (t < 0.601)) || 
                ((t > 0.799) && (t < 0.801)) then
                -200.
            else
                0.
    in
    ((linear @+. non_linear) @+. control ) @+. (p &*. b)
;;

(* Начальное состояние модели 
 * TODO: вынести в конфиг *)
let x_0 = Array2.of_array float64 c_layout [|
    [| 3. +. (atan (-5.)) |];
    [| 4. |];
    [| 4. |]|]
;;


(* Матрица Q функционала *)
let q_r = Array2.of_array float64 c_layout [|
    [| 1.2 ; 0. ; 0. |];
    [| 0. ; 0. ; 0. |];
    [| 0. ; 0. ; 0. |] |]
;;

(* Матрица P функционала *)
let p_r = Array2.of_array float64 c_layout [|
    [| 1. ; 0. ; 0. |];
    [| 0. ; 0. ; 0. |];
    [| 0. ; 0. ; 0. |] |]
;;

(* Целевая функция *)
let z t = 
    Array2.of_array float64 c_layout [|
    [| 3. +. (atan (t *. 10. -. 5.)) |];
    [| 0. |];
    [| 0. |] |]
;;


(* Функция печати результата *)
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

(* Вычисление первого приближения управления*)
let u_func b r h c k_t = 
    let (t, k) = k_t in
    let u'_1_0 = k @*. c in
    let u'_1 = h @-. u'_1_0 in
    let u' = (1./.r) *. (b ><. u'_1) in
    (t, u')
;;

(* Вычисление приближения управления на очередном шаге *)
let u_func_n b r h x k = 
    let u'_0 = k @*. x in
    let u'_1 = h @-. u'_0 in
    let u' = (1./.r) *. (b ><. u'_1) in
    u'
;;

(* Дифференциальное уравнение для расчета значения h *)
let h_diff a b q r f z_f k x t h = 
    let b' = transpose b in
    let h_diff_0 = (1. /. r) &*. ((b @*. b') @*. k) in
    let h_diff_part_0 = a @-. h_diff_0 in 
    let h_diff_part_0' = (transpose h_diff_part_0) @*. h
    and h_diff_part_1 = q @*. (z_f t)
    and h_diff_part_2 = k @*. (f x) in
    (h_diff_part_2 @-. h_diff_part_0') @-. h_diff_part_1
;;

(* Вычисление первого приближения модели *)
let x_new h diff_funct prev u = 
    let (_, _, x)::_ = prev
    and (t, u') = u in
    let x_new = mx_rk4_step (diff_funct u') x t h in
    (t +. h, u', x_new)::prev
;;

(* Вычисление приближения ф-ии h на очередном шаге *)
let h_new step diff_funct prev k_r x_p = 
    let (t, k) = k_r
    and (_, _, x) = x_p
    and (_, h)::_ = prev in
    let h_new =  mx_rk4_step (diff_funct k x) h t step in
    (t +. step, h_new)::prev
;;

(* Вычисление приближения модели на очередном шаге *)
let x_new_n step u_f diff_func prev h_n k_n =
    let (t, u, x)::_  = prev 
    and (_, h_n') = h_n
    and (_, k_n') = k_n in
    let x_new = mx_rk4_step (diff_func u) x t step in
    let u_new = u_f h_n' x_new k_n' in
    (t +. step, u_new, x_new)::prev
;;

(* Вычисление следующего приближения решения *)
let x_next calc_step h_diff_f h_new_f u_next_f 
    x_next_f x_diff_f z_f q b r
    krev k k0 k1 x0 x_prev = 
    (* ВЫчисляем следующее приближение ф-ии h(t)
     * h_new = h_n(0)..h_n(1) *)
    let _::h_new = 
        List.fold_left2 (h_new_f (-.calc_step) h_diff_f)
                ((1., q @*. (z_f 1.))::[])
                ((1., k1)::krev) x_prev in 
    let (_, h_new_0)::_ = h_new in
    (* ВЫчисляем управление при t=0 *)
    let u_new_0 = u_next_f b r h_new_0 x0 k0 in
    (* ВЫчисляем следующее приближение решения x(t)
     * x_n = x_n(1) .. x_n(0) *)
    let x_n =
        List.fold_left2 (x_next_f calc_step (u_next_f b r) x_diff_f)
                ((0., u_new_0, x0)::[])
                h_new k in
    x_n
;;

(* Функция исчисления расстояния между двумя приближениями решения *)
let metric_func a b = 
    let diff acc x y = 
        let (_,_,x') = x
        and (_,_,y') = y in
        let diff =  x' @-. y' in
        let s = diff ><. diff in
        acc +. s
    in
    List.fold_left2 diff 0. a b
;;

(* Реализация метода последовательных приближений *)
let fx eps calc_step metric_f x_func x_0 = 
    (*Вспомогательная ф-ия для реализации рекурсии*)
    let rec fx' i metric_prev x_p =
        (* Вычисляем следующее приближение*)
        let _::x = x_func x_p in
        (* И расстояние от предыдущего приближения *)
        let metric  = (metric_f x_p x) in
        (* И если разница метрик меньше чем эпсилон*)
        match (metric_prev -. metric) < eps with
            (* Возвращаем значение *)
          true  -> Printf.printf "%d\n" i; x
            (* Уходим на второй круг *)
        | false -> fx' (i+1) metric x
    in
    (* Первое приближение *)
    let _::x_1 = x_func x_0 in
    let metric_0 = (metric_f x_0 x_1) in
    fx' 0 metric_0 x_1
;;

(* Главная ф-ия *)
let _ =
    (* Читаем параметры из файла *)
    let in_file = Pervasives.open_in "in.txt" in
    let step = Scanf.fscanf in_file "%f " (fun x -> x) 
    and eps  = Scanf.fscanf in_file "%f " (fun x -> x) 
    and a = Scanf.fscanf in_file "%f " (fun x -> x) 
    and b = Scanf.fscanf in_file "%f " (fun x -> x)
    and d = Scanf.fscanf in_file "%f " (fun x -> x) 
    and sigma = Scanf.fscanf in_file "%f " (fun x -> x)
    and delta = Scanf.fscanf in_file "%f " (fun x -> x) 
    and r_r = Scanf.fscanf in_file "%f " (fun x -> x)
    and _ = q_r.{0,0} <- Scanf.fscanf in_file "%f " (fun x -> x)
    and _ = p_r.{0,0} <- Scanf.fscanf in_file "%f " (fun x -> x) in
    (* Зададим модель *)
    let k_1 = p_r in
    let a_r = a_magnit a d
    and b_r = b_magnit
    and f_r = f_magnit b sigma delta in
    let x_magnit_d = x_magnit_diff a_r b_r f_r 
    and h_magnit_d = h_diff a_r b_r q_r r_r f_r z in
    (* Решим уравнение Риккати
     * k_r = k(0) :: ... :: k(1)  *)
    let k_r = mx_rk4_down (k_magnit_diff a_r b_r q_r r_r) k_1 1. step 0. in
    (* И переврнем решение -- пригодится
     * k_rev = k(1) :: .. :: k(0) *)
    let _::k_rev = List.rev k_r in
    let (_,k_0)::_ = k_r in
    (* Вычисляем нулевое приближение решения*)
    let h_0 = p_r @*. z(1.) in
    (* u_0 = u_0(0)...u_0(1) *)
    let u_0 = List.rev_map (u_func b_r r_r h_0 x_0) k_rev in
    let (_, u0)::_ = u_0 in
    let x_0'= List.fold_left (x_new step x_magnit_d) ((0. , u0, x_0)::[]) u_0 in
    (* Решаем задачу методом последовательных приближений *)
    let x_21' = fx eps step metric_func 
                (x_next step h_magnit_d h_new u_func_n 
                x_new_n x_magnit_d z q_r b_r r_r 
                k_rev k_r k_0 k_1 x_0) x_0' 
    in
    (* И напечатаем результат *)
    List.map (print_u_x) x_21';
;;
