;; Fix colors (like ls --color, etc)
;;; Shell mode
(setq ansi-color-names-vector ; better contrast colors
      ["black" "red4" "green4" "yellow4"
       "blue3" "magenta4" "cyan4" "white"])
(add-hook 'shell-mode-hook 'ansi-color-for-comint-mode-on)

;;Make the prompt read only
(setq comint-prompt-read-only t)

;; turn on dirtrack-mode instead of the shel-dirtrack-mode(default)
(add-hook 'shell-mode-hook
          (lambda ()
            (shell-dirtrack-mode -1)
            (insert "export PS1=\"wanding:\\w$ \"")
            (comint-send-input)
            (dirtrack-mode 1)
	    ))

;; the following list tells the dirtrack mode how to parse the prompt
;; and where to look for the name of the directory
;; the first item is the regular expression of the prompt
;; the secnd item is the group number of the working directory address
(setq-default dirtrack-list '("^wanding:\\([^$]*\\)\\$" 1))

;; make the shell buffer name the directory name
(add-hook 'dirtrack-directory-change-hook
          (lambda ()
            (let ((base-buffer-name (concat "shell-" default-directory "-shell"))
                  (i 1)
                  (full-buffer-name base-buffer-name))
              (while (get-buffer full-buffer-name)
                (setq i (1+ i))
                (setq full-buffer-name (concat base-buffer-name "<" (number-to-string i) ">")))
              (rename-buffer full-buffer-name))))

;; avoid problem in ssh
;; (add-hook 'ssh-mode-hook (lambda () (setq dirtrackp nil)))