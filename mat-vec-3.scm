;;;
;;; Generic 3 vectors and 3X3 matrix operations
;;;
;;; Copyright 2013 Benjamin Silbaugh
;;;
;;; See LICENSE file for redistribution and modification permissions.

(functor (mat-vec-3 (M (add sub mul div sqroot mx)))
  (make-vec
   make-mat
   mag
   vec*sca
   sca*vec
   vec/sca
   vec+vec
   vec-vec
   dot
   cross
   mat*sca
   sca*mat
   mat/sca
   mat*vec
   mat+mat
   mat-mat
   mat*mat)

  (import scheme chicken M)

;;; === Misc Utilities ===

(define (sq x)
  (mul x x))

;;; === Vector and Matrix Data Structures ===

;; Vector in R3
;;
;; The vector can be thought of as a 3-tuple; e.g. (u, v, w). When
;; performing operations with matricies, it will be implicitly identified
;; by the operation as either a 3X1 or 1X3 matrix.
(define-record vec
  elem-1
  elem-2
  elem-3)

;; Matrix in R3X3
(define-record mat
  elem-11 elem-12 elem-13
  elem-21 elem-22 elem-23
  elem-31 elem-32 elem-33)

;;; === Map and Apply Operations ===

;; Evaluates function f with elements of vector v as arguments
(define (vec-apply f v)
  (f (vec-elem-1 v) (vec-elem-2 v) (vec-elem-3 v)))

;; Evaluates function f with elements of matrix a as arguments
(define (mat-apply f a)
  (f (mat-elem-11 a) (mat-elem-12 a) (mat-elem-13 a)
	 (mat-elem-21 a) (mat-elem-22 a) (mat-elem-23 a)
	 (mat-elem-31 a) (mat-elem-32 a) (mat-elem-33 a)))

;; Maps unitary function over vector elements
(define (vec-uni-map f v)
  (make-vec (f (vec-elem-1 v)) (f (vec-elem-2 v)) (f (vec-elem-3 v))))

;; Maps unitary function over matrix elements
(define (mat-uni-map f a)
  (make-mat (f (mat-elem-11 a)) (f (mat-elem-12 a)) (f (mat-elem-13 a))
			(f (mat-elem-21 a)) (f (mat-elem-22 a)) (f (mat-elem-23 a))
			(f (mat-elem-31 a)) (f (mat-elem-32 a)) (f (mat-elem-33 a))))

;; Maps binary function over vector elements
(define (vec-bin-map f u v)
  (make-vec (f (vec-elem-1 u) (vec-elem-1 v))
			(f (vec-elem-2 u) (vec-elem-2 v))
			(f (vec-elem-3 u) (vec-elem-3 v))))

;; Maps binary function over matrix elements
(define (mat-bin-map f a b)
  (make-mat ; -- 1st row ---
            (f (mat-elem-11 a) (mat-elem-11 b))
			(f (mat-elem-12 a) (mat-elem-12 b))
			(f (mat-elem-13 a) (mat-elem-13 b))
			;-- 2nd row --
            (f (mat-elem-21 a) (mat-elem-21 b))
			(f (mat-elem-22 a) (mat-elem-22 b))
			(f (mat-elem-23 a) (mat-elem-23 b))
			;-- 3rd row --
            (f (mat-elem-31 a) (mat-elem-31 b))
			(f (mat-elem-32 a) (mat-elem-32 b))
			(f (mat-elem-33 a) (mat-elem-33 b))))

;;; === Vector Operations ===

;; Computes the magnitude of vector v (euclidean norm)
(define (mag v)
  (let ((u (vec-elem-1 v)) (v (vec-elem-2 v)) (w (vec-elem-3 v)))
	(let ((r (mx u (mx v w))))
	  (mul r (sqroot (add (sq (div u r)) (sq (div v r)) (sq (div w r))))))))

;;; === Vector X Scalar Operations ===

(define (vec*sca v c)
  (vec-uni-map (lambda (x) (mul c x)) v))

(define (sca*vec c v) (vec*sca v c))

(define (vec/sca v c)
  (vec-uni-map (lambda (x) (div x c)) v))

;;; === Vector X Vector Operations ===

(define (vec+vec u v)
  (vec-bin-map add u v))

(define (vec-vec u v)
  (vec-bin-map sub u v))

(define (dot u v)
  (let ((u1 (vec-elem-1 u)) (u2 (vec-elem-2 u)) (u3 (vec-elem-3 u))
		(v1 (vec-elem-1 v)) (v2 (vec-elem-2 v)) (v3 (vec-elem-3 v)))
	(add (mul u1 v1) (mul u2 v2) (mul u3 v3))))

(define (cross u v)
  (let ((u1 (vec-elem-1 u)) (u2 (vec-elem-2 u)) (u3 (vec-elem-3 u))
		(v1 (vec-elem-1 v)) (v2 (vec-elem-2 v)) (v3 (vec-elem-3 v)))
	(make-vec (sub (mul u2 v3) (mul u3 v2))
			  (sub (mul u3 v1) (mul u1 v3))
			  (sub (mul u1 v2) (mul u2 v1)))))

;;; === Matrix Row and Column Operations ===

(define (mat-row-1 a)
  (mat-apply (lambda (a11 a12 a13 
					  a21 a22 a23 
					  a31 a32 a33) (make-vec a11 a12 a13)) a))

(define (mat-row-2 a)
  (mat-apply (lambda (a11 a12 a13 
					  a21 a22 a23 
					  a31 a32 a33) (make-vec a21 a22 a23)) a))

(define (mat-row-3 a)
  (mat-apply (lambda (a11 a12 a13 
					  a21 a22 a23 
					  a31 a32 a33) (make-vec a31 a32 a33)) a))

(define (mat-col-1 a)
  (mat-apply (lambda (a11 a12 a13 
					  a21 a22 a23 
					  a31 a32 a33) (make-vec a11 a21 a31)) a))

(define (mat-col-2 a)
  (mat-apply (lambda (a11 a12 a13 
					  a21 a22 a23 
					  a31 a32 a33) (make-vec a12 a22 a32)) a))

(define (mat-col-3 a)
  (mat-apply (lambda (a11 a12 a13 
					  a21 a22 a23 
					  a31 a32 a33) (make-vec a13 a23 a33)) a))

;;; === Matrix X Scalar Operations ===

(define (mat*sca a c)
  (mat-bin-map (lambda (aij) (mul c aij)) a))

(define (sca*mat c a) (mat*sca a c))

(define (mat/sca a c)
  (mat-bin-map (lambda (aij) (div aij c)) a))

;;; === Matrix X Vector Operations ===

;; Matrix-vector multiplication
;;
;; Computes [a][v], where [a] and [v] are 3X3 and 3X1 matricies,
;; respectively.
(define (mat*vec a u)
  (make-vec (dot (mat-row-1 a) u) 
			(dot (mat-row-2 a) u) 
			(dot (mat-row-3 a) u)))

;; Vector-matrix multiplication
;;
;; Computes [v][a], where [v] and [a] are 1X3 and 3X3 matricies,
;; respectively.
(define (vec*mat u a)
  (make-vec (dot u (mat-col-1 a)) 
			(dot u (mat-col-2 a)) 
			(dot u (mat-col-3 a))))

;; === Matrix X Matrix Operations ===

(define (mat+mat a b)
  (mat-bin-map add a b))

(define (mat-mat a b)
  (mat-bin-map sub a b))

(define (mat*mat a b)
  (let ((a1 (mat-row-1 a)) (a2 (mat-row-2 a)) (a3 (mat-row-3 a))
		(b1 (mat-col-1 b)) (b2 (mat-col-2 b)) (b3 (mat-col-3 b)))
	(make-mat (dot a1 b1) (dot a1 b2) (dot a1 b3)
			  (dot a2 b1) (dot a2 b2) (dot a2 b3)
			  (dot a3 b1) (dot a3 b2) (dot a3 b3))))

); end functor
