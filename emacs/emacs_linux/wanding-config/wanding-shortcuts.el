;; f1 is help
;; f2 is menu
;; f3-f4 is macro
(global-set-key (kbd "<f5>") 'beginning-of-buffer)
(global-set-key (kbd "<f6>") 'cua-paste)
(global-set-key (kbd "<f7>") 'goto-line)
(global-set-key (kbd "<f8>") 'other-window)
;; (global-set-key (kbd "<f9>") 'find-file)
;; (global-set-key (kbd "<f9>") 'delete-other-windows)
(global-set-key (kbd "<f9>") 'my-switch-to-other-buffer)
(global-set-key (kbd "<f10>") 'helm-recentf)
;; (global-set-key (kbd "<f11>") 'toggle-fullscreen)
;; (global-set-key (kbd "<f11>") 'switch-to-buffer)
(global-set-key (kbd "<f11>") 'dabbrev-expand)
(global-set-key (kbd "<f12>") 'save-buffer)

(global-set-key (kbd "M-a") 'unfill-paragraph)

;; (global-set-key (kbd "C-x f") 'toggle-fullscreen)
(global-set-key [XF86Forward] 'other-window)
(global-set-key (kbd "M-<up>") 'move-line-upward)
(global-set-key (kbd "M-<down>") 'move-line-downward)
;; this setup meet the situation where M-v is scroll-down
(global-set-key (kbd "M-c") 'scroll-up)
(global-set-key (kbd "M-i") 'kill-whole-line)
;; C-x C-b invoke buffer list in the current window (buffer-menu instead of list-buffer)
(global-set-key "\C-x\C-b" 'buffer-menu)

(global-set-key (kbd "<pause>") 'dabbrev-expand)
(global-set-key (kbd "<insert>") 'other-window)
(global-set-key (kbd "<Scroll_Lock>") 'set-mark-command)

;;;;; Make page up and page down a whole lot nicer
;; (global-set-key "\C-v"	   'pager-page-down)
;; (global-set-key [next] 	   'pager-page-down)
;; (global-set-key "\ev"	   'pager-page-up)
;; (global-set-key [prior]	   'pager-page-up)
;; (global-set-key '[M-up]    'pager-row-up)
;; (global-set-key '[M-kp-8]  'pager-row-up)
;; (global-set-key '[M-down]  'pager-row-down)
;; (global-set-key '[M-kp-2]  'pager-row-down)
