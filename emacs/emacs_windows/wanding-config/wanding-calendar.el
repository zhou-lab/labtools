;; Fix foolish calendar-mode scrolling.
(add-hook 'calendar-load-hook
	  '(lambda ()
	     (setq mark-holidays-in-calendar t)
	     (define-key calendar-mode-map ">" 'scroll-calendar-left)
	     (define-key calendar-mode-map "<" 'scroll-calendar-right)
	     (define-key calendar-mode-map "\C-x>" 'scroll-calendar-left)
	     (define-key calendar-mode-map "\C-x<" 'scroll-calendar-right)))
