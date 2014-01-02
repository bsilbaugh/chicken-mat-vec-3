;;;
;;; Test suite for mat-vec-3
;;;
;;; The test suite mostly consists of "sanity checks" (i.e. verifying basic
;;; algebraic identities hold), plus a few sample calculations. That said,
;;; the test suite is not perfect: some programming errors may still manage
;;; to remain undetected.
;;;

(use test)

(use mat-vec-3)

(module fp-mat-vec-3 = mat-vec-3
  (import scheme chicken)
  (define add fp+)
  (define sub fp-)
  (define mul fp*)
  (define div fp/)
  (define sqroot fpsqrt)
  (define (maxabs x y)
	(fpmax (fpabs x) (fpabs y))))

(import fp-mat-vec-3)

;; === Test values ===

(define zero-vec (make-vec 0.0 0.0 0.0))
(define u (make-vec 0.1 0.2 0.3))
(define v (make-vec 0.01 0.02 0.03))
(define w (make-vec 0.001 -0.02 0.003))

(define zero-mat (make-mat 0.0 0.0 0.0
						   0.0 0.0 0.0
						   0.0 0.0 0.0))

(define id3 (make-mat 1.0 0.0 0.0
					  0.0 1.0 0.0
					  0.0 0.0 1.0))

(define a (make-mat  0.102  0.143  0.245
					 0.842 -0.315 -0.143
					-0.452  0.644  0.952))

(define b (make-mat 0.432 0.743 0.836
					0.472 0.844 0.824
					0.457 0.379 0.867))

;; === Helper utitlities ===

(define (close? eps x y)
  (> eps (abs (- x y))))

(define (vec-close? eps u v)
  (and (close? eps (vec-elem-1 u) (vec-elem-1 v))
	   (close? eps (vec-elem-2 u) (vec-elem-2 v))
	   (close? eps (vec-elem-3 u) (vec-elem-3 v))))

(define (mat-close? eps a b)
  (and (close? eps (mat-elem-11 a) (mat-elem-11 b))
	   (close? eps (mat-elem-12 a) (mat-elem-12 b))
	   (close? eps (mat-elem-13 a) (mat-elem-13 b))
	   (close? eps (mat-elem-21 a) (mat-elem-21 b))
	   (close? eps (mat-elem-22 a) (mat-elem-22 b))
	   (close? eps (mat-elem-23 a) (mat-elem-23 b))
	   (close? eps (mat-elem-31 a) (mat-elem-31 b))
	   (close? eps (mat-elem-32 a) (mat-elem-32 b))
	   (close? eps (mat-elem-33 a) (mat-elem-33 b))))

;; === Vector X scalar operations ===

(test-group "mag"
  (test 0.0 (mag (make-vec 0.0 0.0 0.0)))
  (test 0.1 (mag (make-vec 0.1 0.0 0.0)))
  (test 0.1 (mag (make-vec 0.0 0.1 0.0)))
  (test 0.1 (mag (make-vec 0.0 0.0 0.1)))
  (test 0.1 (mag (make-vec -0.1  0.0  0.0)))
  (test 0.1 (mag (make-vec  0.0 -0.1  0.0)))
  (test 0.1 (mag (make-vec  0.0  0.0 -0.1))))

(test-group "vec*sca"
  (test-assert (equal? v (vec*sca v 1.0)))
  (test-assert (equal? (make-vec 0.0 0.0 0.0) (vec*sca v 0.0)))
  (test-assert (equal? (make-vec 0.05 0.01 0.005) 
					   (vec*sca (make-vec 0.1 0.02 0.01) 0.5))))

(test-group "sca*vec"
  (test-assert (equal? (make-vec 0.05 0.01 0.005) 
					   (sca*vec 0.5 (make-vec 0.1 0.02 0.01)))))

(test-group "vec/sca"
  (test-assert (equal? (make-vec 0.05 0.01 0.005) 
					   (vec/sca (make-vec 0.1 0.02 0.01) 2.0))))


;; === Vector X vector operations ===

(test-group "vec+vec"
  (test-assert (equal? v (vec+vec zero-vec v)))
  (test-assert (equal? v (vec+vec v zero-vec)))
  (test-assert (not (equal? zero-vec (vec+vec u v))))
  (test-assert (vec-close? 1.0E-6 
			               (make-vec 0.1 -0.1 0.2) 
					       (vec+vec (make-vec 0.05  0.2  0.4)
								    (make-vec 0.05 -0.3 -0.2)))))

(test-group "vec-vec"
  (test-assert (equal? v (vec-vec v zero-vec)))
  (test-assert (not (equal? zero-vec (vec-vec u v))))
  (test-assert (vec-close? 1.0E-6 
			               (make-vec 0.0 0.5 0.6) 
					       (vec-vec (make-vec 0.05  0.2  0.4)
								    (make-vec 0.05 -0.3 -0.2)))))

(test-group "dot"
  (test 1.0 (dot (make-vec 1.0 0.0 0.0) (make-vec 1.0 0.0 0.0)))
  (test 0.0 (dot (make-vec 1.0 0.0 0.0) (make-vec 0.0 1.0 0.0)))
  (test 0.0 (dot (make-vec 1.0 0.0 0.0) (make-vec 0.0 0.0 1.0)))
  (test 0.0 (dot (make-vec 0.0 1.0 0.0) (make-vec 1.0 0.0 0.0)))
  (test 1.0 (dot (make-vec 0.0 1.0 0.0) (make-vec 0.0 1.0 0.0)))
  (test 0.0 (dot (make-vec 0.0 1.0 0.0) (make-vec 0.0 0.0 1.0)))
  (test 0.0 (dot (make-vec 0.0 0.0 1.0) (make-vec 1.0 0.0 0.0)))
  (test 0.0 (dot (make-vec 0.0 0.0 1.0) (make-vec 0.0 1.0 0.0)))
  (test 1.0 (dot (make-vec 0.0 0.0 1.0) (make-vec 0.0 0.0 1.0))))

(test-group "cross"
  (test-assert (equal? (make-vec 0.0 0.0 0.0)
					   (cross (make-vec 1.0 0.0 0.0) 
							  (make-vec 1.0 0.0 0.0))))
  (test-assert (equal? (make-vec 0.0 0.0 1.0)
					   (cross (make-vec 1.0 0.0 0.0) 
							  (make-vec 0.0 1.0 0.0))))
  (test-assert (equal? (make-vec 0.0 -1.0 0.0)
					   (cross (make-vec 1.0 0.0 0.0) 
							  (make-vec 0.0 0.0 1.0))))
  (test-assert (equal? (make-vec 0.0 0.0 -1.0)
					   (cross (make-vec 0.0 1.0 0.0) 
							  (make-vec 1.0 0.0 0.0))))
  (test-assert (equal? (make-vec 0.0 0.0 0.0)
					   (cross (make-vec 0.0 1.0 0.0) 
							  (make-vec 0.0 1.0 0.0))))
  (test-assert (equal? (make-vec 1.0 0.0 0.0)
					   (cross (make-vec 0.0 1.0 0.0) 
							  (make-vec 0.0 0.0 1.0))))
  (test-assert (equal? (make-vec 0.0 1.0 0.0)
					   (cross (make-vec 0.0 0.0 1.0) 
							  (make-vec 1.0 0.0 0.0))))
  (test-assert (equal? (make-vec -1.0 0.0 0.0)
					   (cross (make-vec 0.0 0.0 1.0) 
							  (make-vec 0.0 1.0 0.0))))
  (test-assert (equal? (make-vec 0.0 0.0 0.0)
					   (cross (make-vec 0.0 0.0 1.0) 
							  (make-vec 0.0 0.0 1.0)))))

;; === Matrix X scalar operations ===

(test-group "mat*sca"
  (test-assert (equal? zero-mat (mat*sca a 0.0)))
  (test-assert (equal? a (mat*sca a 1.0)))
  (test-assert (mat-close? 1.0E-6 
						   (make-mat 0.11 0.12 0.13
									 0.21 0.22 0.23
									 0.31 0.32 0.33)
						   (mat*sca
							(make-mat 1.1 1.2 1.3
									  2.1 2.2 2.3
									  3.1 3.2 3.3) 0.1))))

(test-group "sca*mat"
  (test-assert (equal? (mat*sca a 0.1) (sca*mat 0.1 a))))

(test-group "mat/sca"
  (test-assert (equal? a (mat/sca a 1.0)))
  (test-assert (mat-close? 1.0E-6 
						   (make-mat 0.11 0.12 0.13
									 0.21 0.22 0.23
									 0.31 0.32 0.33)
						   (mat/sca
							(make-mat 1.1 1.2 1.3
									  2.1 2.2 2.3
									  3.1 3.2 3.3) 10.0))))

;; === Matrix X vector operations ===

(test-group "mat*vec"
  (test-assert (equal? zero-vec (mat*vec a zero-vec)))
  (test-assert (equal? zero-vec (mat*vec zero-mat u)))
  (test-assert (equal? u (mat*vec id3 u)))
  (test-assert (equal? (mat-col-1 a) (mat*vec a (make-vec 1.0 0.0 0.0))))
  (test-assert (equal? (mat-col-2 a) (mat*vec a (make-vec 0.0 1.0 0.0))))
  (test-assert (equal? (mat-col-3 a) (mat*vec a (make-vec 0.0 0.0 1.0)))))

(test-group "mat+mat"
  (test-assert (equal? a (mat+mat a zero-mat)))
  (test-assert (mat-close? 1.0E-6
						   (make-mat 11.11 12.12 13.13
									 21.21 22.22 23.23
									 31.31 32.32 33.33)
						   (mat+mat 
							(make-mat 11.0 12.0 13.0
									  21.0 22.0 23.0
									  31.0 32.0 33.0)
							(make-mat 0.11 0.12 0.13
									  0.21 0.22 0.23
									  0.31 0.32 0.33)))))

(test-group "mat-mat"
  (test-assert (equal? a (mat-mat a zero-mat)))
  (test-assert (mat-close? 1.0E-6 
						   (mat+mat a (mat*sca b -1.0))
						   (mat-mat a b))))

(test-group "mat*mat"
  (test-assert (equal? zero-mat (mat*mat a zero-mat)))
  (test-assert (equal? a (mat*mat a id3))))
						   

(test-exit)
