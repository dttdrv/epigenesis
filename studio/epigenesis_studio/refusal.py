class Refusal(ValueError):
    """A refusal the page can phrase in the viewer's language: `code` and `params` travel with the English text."""

    def __init__(self, code, text, **params):
        super().__init__(text)
        self.code, self.params = code, params
