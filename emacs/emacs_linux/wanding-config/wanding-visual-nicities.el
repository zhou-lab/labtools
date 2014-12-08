(add-to-list 'custom-theme-load-path "~/.emacs.d/wanding-config/themes")
;; see http://orgmode.org/worg/org-color-themes.html for a list
;; (load-theme 'leuven t)
;; (load-theme 'cyberpunk t)
;; (load-theme 'zenburn t)
;; (load-theme 'wombat t)
;; (load-theme 'ample t)
;; (load-theme 'manoj-dark t)
;; (load-theme 'tangotango t)
;; (load-theme 'sanityinc-solarized-light t)
;; (load-theme 'misterioso t)
;; (load-theme 'dichromacy t)

(setq default-frame-alist
      '(
	(width . 110)
	(height . 35)
	)
      )

(transient-mark-mode 1)
;; (require 'color-theme)
; (load-library "wanding-color-theme")
(global-font-lock-mode 1)
; (color-theme-initialize)
; (color-theme-ryan)
;; (setq tab-width 2)
;; (setq-default tab-width 2)
;; (setq indent-tabs-mode nil)
;; (setq c-basic-offset 2)
;; (setq perl-indent-level 2)
;; (setq bash-indent-level 2)


;; put tabbar
;; (require 'tabbar)
(tabbar-mode t)

;; (menu-bar-mode -1)
;; (tool-bar-mode -1)
;; (scroll-bar-mode -1)


;; display line and column number
(setq column-number-mode t)
(setq line-number-mode t)

;; do not display startup message
(setq inhibit-startup-message t)
    
;; disable beep
(setq visible-bell t)

;; show matched parentheses
(show-paren-mode t)

;; the frame title has only the filename
(setq frame-title-format
(concat "%b"))


;;;;; font settings ;;;;;

;; (set-default-font "MONACO 8")
;; (set-default-font "Yahei Consolas Hybrid 10")
;; (set-fontset-font "fontset-default" 'gb18030 '("Yahei Consolas Hybrid 11" . "unicode-bmp"))
;; (set-default-font "WenQuanYi Zen Hei Mono 10")

;; If you open a new frame, and the font is not the specified one,
;; use the following instead of set-default-font
;; (add-to-list 'default-frame-alist '(font . "Yahei Consolas Hybrid 9"))
;; (add-to-list 'default-frame-alist '(font . "Monaco-10"))
;; (add-to-list 'default-frame-alist '(font . "Monaco-8"))
;; (add-to-list 'default-frame-alist '(font . "Lucida Grande 10"))
;; (add-to-list 'default-frame-alist '(font . "Lucida Grande 10"))
;; (add-to-list 'default-frame-alist '(font . "Nimbus Mono L 10"))
;; (add-to-list 'default-frame-alist '(font . "FreeMono 13"))


;; (add-to-list 'default-frame-alist '(font . "FZZY - Kelvin 11"))

;; (add-to-list 'default-frame-alist '(font . "WenQuanYi Micro Hei Mono 9"))
;; (add-to-list 'default-frame-alist '(font . "Monospace 13"))
(add-to-list 'default-frame-alist '(font . "Droid Sans Mono 13"))
;; (add-to-list 'default-frame-alist '(font . "WenQuanYi Micro Hei 9"))


;; set mouse scroll
(setq mouse-wheel-scroll-amount '(2))
(setq mouse-wheel-progressive-speed nil)

(defun toggle-line-spacing ()
  "Toggle line spacing between 1 and 5 pixels."
  (interactive)
  (if (eq line-spacing 1)
      (setq-default line-spacing 0)
    (setq-default line-spacing 1))
  )


(defun toggle-fullscreen (&optional f)
  (interactive)
  (let ((current-value (frame-parameter nil 'fullscreen)))
    (set-frame-parameter nil 'fullscreen
			 (if (equal 'fullboth current-value)
			     (if (boundp 'old-fullscreen) old-fullscreen nil)
			   (progn (setq old-fullscreen current-value)
				  'fullboth)))))


;Show column numbers
(column-number-mode 1)
(setq-default fill-column 72)
(setq auto-fill-mode 1)
;Show what's being selected
(transient-mark-mode 1)
;Line by line scrolling
(setq scroll-step 1)
(setq inhibit-startup-message t)

;; (setq color-theme-load-all-themes nil)

;; (require 'color-theme-tangotango)

;; select theme - first list element is for windowing system, second is for console/terminal
;; Source : http://www.emacswiki.org/emacs/ColorTheme#toc9
;; (require 'color-theme)
;; (color-theme-initialize)
;; (color-theme-jsc-light2)
;; (color-theme-clarity)
;; (color-theme-pok-wog)
;; (color-theme-andreas)
;; (color-theme-arjen)
;; (color-theme-aalto-light)
;; (color-theme-comidia)
;; (color-theme-infodoc)

;; ;; default-start
;; (funcall (lambda (cols)
;;     	   (let ((color-theme-is-global nil))
;;     	     (eval 
;;     	      (append '(if (window-system))
;;     		      (mapcar (lambda (x) (cons x nil)) 
;;     			      cols)))))
;;     	 color-theme-choices)

;; ;; test for each additional frame or console
;; (require 'cl)
;; (fset 'test-win-sys 
;;       (funcall (lambda (cols)
;;     		 (lexical-let ((cols cols))
;;     		   (lambda (frame)
;;     		     (let ((color-theme-is-global nil))
;; 		       ;; must be current for local ctheme
;; 		       (select-frame frame)
;; 		       ;; test winsystem
;; 		       (eval 
;; 			(append '(if (window-system frame)) 
;; 				(mapcar (lambda (x) (cons x nil)) 
;; 					cols)))))))
;;     	       color-theme-choices ))
;; ;; hook on after-make-frame-functions
;; (add-hook 'after-make-frame-functions 'test-win-sys)

;; ;;(color-theme-tangotango)
;; -- appearance --

;; (tool-bar-mode -1)
;; (menu-bar-mode -1)


; turn on syntax coloring just in case it hasn't been turned on by default.
;; (global-font-lock-mode 1)

; display line numbers in margin (fringe)
;; (global-linum-mode 1)


; set background color
;; (setq default-frame-alist
;;      '( (cursor-color . "white" )
;;         (background-color . "black")
;;         (foreground-color . "white")
;;         )
;;      )


;; (defun set-frame-size-according-to-resolution ()
;;   (interactive)
;;   (if window-system
;;       (progn
;; 	;; use 120 char wide window for largeish displays
;; 	;; and smaller 80 column windows for smaller displays
;; 	;; pick whatever numbers make sense for you
;; 	;; (if (> (x-display-pixel-width) 1500)
;; 	;;    (add-to-list 'default-frame-alist (cons 'width 120))
;; 	(add-to-list 'default-frame-alist (cons 'width 80))
;; 	;; for the height, subtract a couple hundred pixels
;; 	;; from the screen height (for panels, menubars and
;; 	;; whatnot), then divide by the height of a char to
;; 	;; get the height we want
;; 	(add-to-list 'default-frame-alist 
;; 		     (cons 'height (/ (- (x-display-pixel-height) 400) (frame-char-height)))))))
;;(set-frame-size-according-to-resolution)
;; (desktop-save-mode t)
