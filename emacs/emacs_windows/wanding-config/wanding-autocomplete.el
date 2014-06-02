(require 'auto-complete)
(global-auto-complete-mode t)

;; use tab to complete
(define-key ac-complete-mode-map "\t" 'ac-complete)
(define-key ac-complete-mode-map "\r" nil)

;; start autocomplete after typing 3 characters
(setq ac-auto-start 3)
