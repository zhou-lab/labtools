(autoload 'matlab-mode "matlab" "Enter MATLAB mode." t)
(setq auto-mode-alist (cons '("\\.m\\'" . matlab-mode) auto-mode-alist))
(autoload 'matlab-shell "matlab" "Interactive MATLAB mode." t)

(setq matlab-indent-function-body t)  ; if you want function bodies indented
(setq matlab-verify-on-save-flag nil) ; turn off auto-verify on save
(defun my-matlab-mode-hook ()
  (setq fill-column 76))              ; where auto-fill should wrap
(add-hook 'matlab-mode-hook 'my-matlab-mode-hook)
(defun my-matlab-shell-mode-hook ()
  '())

(add-hook 'matlab-shell-mode-hook 'my-matlab-shell-mode-hook)