;; perl indentation
(fset 'perl-mode 'cperl-mode)
(setq cperl-indent-level 4
      cperl-close-paren-offset -4
      cperl-continued-statement-offset 4
      cperl-indent-parens-as-block t
      cperl-tab-always-indent t)
(add-hook 'cperl-mode-hook '(lambda () (progn (pabbrev-mode)) ))