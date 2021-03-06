;; dired mode setup
(require 'dired-single)
(require 'dired-details)
(require 'dired-details+)
(require 'dired-x)

;; (add-hook 'dired-load-hook '(lambda () (require 'dired-x)))
;; (add-hook 'dired-load-hook
;; 	  (lambda()
;; 	    (load "dired-x")
;; 	    (dired-omit-mode 1)))

(dired-details-install)
;; (dired-details-hide)

(add-hook 'dired-mode-hook
	  (lambda()
	    (setq dired-omit-files-p t)
	    ;; omit the dot-files
	    (setq dired-omit-files (concat dired-omit-files "\\|^\\..+$"))
	    (setq dired-omit-extensions '(".pyc" "~" ".bak"))
	    (dired-omit-mode 1)))

;Make sure each dired buffer doesn't spawn new dired buffers
;; (defun my-dired-init ()
;;   "Bunch of stuff to run for dired, either immediately or when it's
;;   loaded."
;;   ;; <add other stuff here>
;;   (define-key dired-mode-map [return] 'joc-dired-single-buffer)
;;   (define-key dired-mode-map [mouse-1] 'joc-dired-single-buffer-mouse)
;;   (define-key dired-mode-map "^"
;;     (function
;;      (lambda nil (interactive) (joc-dired-single-buffer "..")))))
;; ;; if dired's already loaded, then the keymap will be bound
;; (if (boundp 'dired-mode-map)
;;     ;; we're good to go; just add our bindings
;;     (my-dired-init)
;;   ;; it's not loaded yet, so add our bindings to the load-hook
;;   (add-hook 'dired-load-hook 'my-dired-init))
;; ;;enable recursive deletion of dirs, but doubly ask if it's not empty.
;; (setq dired-recursive-deletes 'top)
