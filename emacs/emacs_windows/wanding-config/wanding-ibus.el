(require 'ibus)
;: dependency:
;; apt-get install python-xlib
(add-hook 'after-init-hook 'ibus-mode-on)
(add-hook 'after-make-frame-functions
	  (lambda (new-frame)
	    (select-frame new-frame)
	    (ibus-mode-on)))
;; use C-SPC for Set Mark
(ibus-define-common-key ?\C-\s nil)
;; use C-/ for Undo
(ibus-define-common-key ?\C-/ nil)
;; Change cursor color depending on IBus status
(setq ibus-cursor-color '("red" "blue" "limegreen"))

(global-set-key (kbd "C-z") 'ibus-toggle)
